## **一、TFHEpp 新版的变化**

1. **接口调整**

   - 原来 sk.key.lvl2 这样的成员没了，现在要用 sk.key.get<lvl2param>()。

     ~~~c++
     #include <tfhe++.hpp>
     using namespace TFHEpp;
     
     SecretKey sk;
     
     // ✅ 批量加密
     std::vector<uint8_t> bits = {0,1,1,0};
     auto cts = bootsSymEncrypt<lvl2param>(bits, sk);
     
     // ✅ 批量解密
     auto dec = bootsSymDecrypt<lvl2param>(cts, sk);
     // dec 是 vector<uint8_t>，与 bits 对应
     ~~~

     

   - bootsSymEncrypt / bootsSymDecrypt 变成了**批量接口**（输入输出 vector<uint8_t>），不再支持原先的一位直接加密/解密。

2. **单比特操作方式**

   - 如果要逐位加密，需要用 tlweSymEncrypt 和 tlweSymDecrypt。
   - 明文 0/1 要手动映射成 -μ / +μ，而不是直接丢个 uint8_t 给它。

   ~~~c++
   #include <tfhe++.hpp>
   using namespace TFHEpp;
   
   template<class P>
   TLWE<P> encrypt_bit(uint8_t b, const SecretKey& sk) {
       TLWE<P> c;
       const auto mu = P::μ;
       const auto m  = b ? +mu : -mu;
       // 注意：key 通过 get<P>() 取
       tlweSymEncrypt<P>(c, m, P::α, sk.key.get<P>());
       return c;
   }
   
   template<class P>
   uint8_t decrypt_bit(const TLWE<P>& c, const SecretKey& sk) {
       return tlweSymDecrypt<P>(c, sk.key.get<P>()); // 返回 0/1
   }
   ~~~

   

3. **GateBootstrapping 模板签名更严格**

   - 不再能随便写 GateBootstrapping(res, tlwe, ek)，需要明确匹配 bkP, iksP 模板参数，或者通过包装好的函数来走。

   



## **二、对应修改**



1. **HomGate 系列**

   - 改了 PBS_Single / PBS_Multi，把 EvalKey& 改成了 const EvalKey&，然后用 const_cast 拿到 ek.get<T>() 的非 const 指针。
   - 这样外部接口保持干净，但内部还是能取出 KeySwitchingKey 和 BootstrappingKeyFFT。
   - HomGate 本身逻辑没变，只是调用方式跟着新版模板接口走。

   ~~~c++
   #pragma once
   #include <mulfft.hpp>
   #include <tfhe++.hpp>
   
   namespace TFHEpp {
   
   // 🔴 关键：inline，避免 ODR 多重定义
   inline void PBS_Single(TLWE<lvl2param> &res,
                          const TLWE<lvl2param> &tlwe,
                          const EvalKey &ek,
                          const std::vector<std::vector<int>> &lut)
   {
       TLWE<lvl0param> tlwe0;
   
       // 对内“去 const”只为了 get<T>() 的返回 shared_ptr 解引用
       auto& ek_nc   = const_cast<EvalKey&>(ek);
       auto  ksk_sp  = ek_nc.get<KeySwitchingKey<lvl20param>>();
       auto& ksk     = *ksk_sp;
       IdentityKeySwitch<lvl20param>(tlwe0, tlwe, ksk);
   
       Polynomial<lvl2param> poly{};
       constexpr int padding = 128;
       for (int j = 0; j < padding; j++) poly[j] = lvl2param::μ;
   
       std::array<typename lvl2param::T, 16> t{};
       const int t_size = (int)lut[0].size();
       for (int i = 0; i < t_size; i++) t[(i + 16 - t_size)] = (lut[0][i]);
   
       for (int i = 0; i < 16; i++)
           for (int j = 0; j < padding; j++)
               poly[i * padding + j] = (t[i] * 2 - 1) * lvl2param::μ;
   
       TRLWE<lvl02param::targetP> acc;
       auto  bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
       auto& bkfft    = *bkfft_sp;
       BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);
   
       SampleExtractIndex<lvl02param::targetP>(res, acc, 0);
   }
   
   // 🔴 关键：inline，避免 ODR 多重定义
   inline void PBS_Multi(std::vector<TLWE<lvl2param>> &res,
                         const TLWE<lvl2param> &tlwe,
                         const EvalKey &ek,
                         const std::vector<std::vector<int>> &lut)
   {
       const int outn   = (int)lut.size();
       const int t_size = (int)lut[0].size();
   
       TLWE<lvl0param> tlwe0;
       auto& ek_nc  = const_cast<EvalKey&>(ek);
       auto  ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
       auto& ksk    = *ksk_sp;
       IdentityKeySwitch<lvl20param>(tlwe0, tlwe, ksk);
   
       Polynomial<lvl2param> poly{};
       constexpr int padding = 128;
       for (int j = 0; j < padding; j++) poly[j] = lvl2param::μ;
   
       TRLWE<lvl02param::targetP> acc, cres;
       auto  bkfft_sp = ek_nc.get<BootstrappingKeyFFT<lvl02param>>();
       auto& bkfft    = *bkfft_sp;
       BlindRotate<lvl02param>(acc, tlwe0, bkfft, poly);
   
       res.resize(outn);
       for (int k = 0; k < outn; k++) {
           std::array<typename lvl2param::T, lvl2param::n> t{};
           for (int i = 0; i < t_size; i++)
               t[(i + 16 - t_size) * padding] = (lut[k][i] * 2 - 1);
           for (int j = 0; j < lvl2param::k + 1; j++)
               PolyMul<lvl2param>(cres[j], acc[j], t);
           SampleExtractIndex<lvl02param::targetP>(res[k], cres, 0);
       }
   }
   
   } // namespace TFHEpp
   ~~~

   

2. **密钥调用**

   - 全部把 sk.key.lvl2 换成了 sk.key.get<lvl2param>()。
   - 这点是最普遍的修改，不然编译直接过不去。

   

3. **测试程序 (test_NOT / test_SUB / test_XOR)**

   - 全部把 bootsSymEncrypt / bootsSymDecrypt 单比特调用换成了 tlweSymEncrypt / tlweSymDecrypt。
   - 映射规则：p ? μ : -μ，噪声用 param::α。
   - 这样每一位都能正确加解密。
   - 还有一个小细节，把测试里变量名 and_cts 统一改成 xor_cts，避免语义混乱。

   



