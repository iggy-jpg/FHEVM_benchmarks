

#include "HomOP.h"
#include "HomUtil.h"
#include "timer.h"
#include <iostream>
#include <random>

using namespace TFHEpp;

int main()
{
    std::cout << "TLWE<lvl2param>=" << sizeof(TLWE<lvl2param>) << "B  "
              << "TLWE<lvl0param>=" << sizeof(TLWE<lvl0param>) << "B\n";
    std::cout << "CAT=16 栈占用: in=" << 16 * sizeof(TLWE<lvl2param>)
              << "B out=" << 16 * sizeof(TLWE<lvl0param>) << "B 合计="
              << (16 * sizeof(TLWE<lvl2param>) + 16 * sizeof(TLWE<lvl0param>)) / 1024
              << "KB\n";

    SecretKey sk;
    EvalKey ek;
    ek.emplaceiksk<lvl20param>(sk);

    
    const int N = 40;
    std::vector<TLWE<lvl2param>> ins(N);
    std::default_random_engine rng(12345);
    for (int i = 0; i < N; ++i) {
        std::vector<TLWE<lvl2param>> bits(4);
        for (int j = 0; j < 4; ++j)
            bits[j] = bootsSymEncrypt<lvl2param>({(uint8_t)(rng() & 1)}, sk)[0];
        Linear<lvl2param>(ins[i], bits, {1, 2, 4, 8}, 16);
    }

    
    std::vector<TLWE<lvl0param>> ref(N);
    {
        auto &ek_nc = const_cast<EvalKey &>(ek);
        auto ksk_sp = ek_nc.get<KeySwitchingKey<lvl20param>>();
        double t = 0;
        { AutoTimer x(&t);
          for (int i = 0; i < N; ++i)
              IdentityKeySwitch<lvl20param>(ref[i], ins[i], *ksk_sp); }
        std::cout << "逐个 IKS: " << t / N << " ms/个\n";
    }

    
    for (int T : {1, 8}) {
        std::vector<TLWE<lvl0param>> got;
        double t = 0;
        { AutoTimer x(&t); BatchKeySwitch20(got, ins, ek, T); }
        size_t diff = 0;
        for (int i = 0; i < N; ++i)
            for (size_t j = 0; j < ref[i].size(); ++j)
                if (ref[i][j] != got[i][j]) ++diff;
        std::cout << "批量 KS (T=" << T << "): " << t / N << " ms/个  位差异字数="
                  << diff << (diff == 0 ? "  [位精确一致 ✓]" : "  [不一致 ✗]") << "\n";
    }
    return 0;
}
