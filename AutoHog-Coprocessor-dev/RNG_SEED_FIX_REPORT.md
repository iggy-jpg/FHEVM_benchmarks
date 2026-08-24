# AutoHog-Coprocessor 随机源与相关安全问题修复报告

本文说明 AutoHog-Coprocessor 中已经完成的随机源、公钥生成、加密随机量、CRS 及外部
输入边界修复。报告以原 `RNG_AUDIT.pdf` 的发现为基准，同时纳入修复过程中发现并处理的
64-bit Torus 噪声缩放、TFHEpp BLAKE3 随机流重叠、序列化、并发错误传播和密文对象
边界问题。

## 1. 修复结果总表

### 1.1 原 `RNG_AUDIT.pdf` 问题对照

| 原编号 | 原报告问题 | 当前状态 | 实际处理 |
|---|---|---|---|
| R1 | 公开种子同时派生公开矩阵和秘密噪声，可从公开 seed 重放噪声 | **已修复** | 公开 seed 只进入 AES-CTR 公开流；秘密误差使用独立 `SecureRng` |
| R2 | seeded RLWE 公钥噪声来自约 31-bit LCG 状态 | **已修复** | `FheKey::keygen()` 的 lvl1/lvl2 公钥生成调用点分别创建 `SecureRng`，外部 `URBG&` 参数实际消费该安全引擎 |
| R3 | ZK 加密的 `r/e0/e1` 来自约 31-bit LCG | **已修复** | ZK 加密入口改用 `SecureRng`，proof seed 统一从安全熵源取满 32 字节 |
| R4 | 非 ZK 普通/V2 加密的 `v/e0/e1` 来自弱引擎 | **已修复** | 普通和 V2 加密入口均改用 `SecureRng` |
| R5 | 短状态 LCG 在约数百次加密后可能出现随机流重叠 | **已修复原根因** | 安全敏感路径不再使用 `minstd_rand/default_random_engine`，每次加密从 CSPRNG 取得新随机流 |
| R6 | 公开 seed 名义为 64 bit，实际只有一次约 32-bit 熵输入 | **已修复** | 公开 seed 改为固定 16 字节，即 128 bit，并一次取满全部字节 |
| R7 | 使用 `mt19937_64` 展开公开 seed，不是密码学伪随机展开 | **已修复** | 改为 AES-128 计数流，使用固定 counter 编码、用途域和互不重叠的流区间 |
| R8 | CRS 陷门 `alpha` 及幂次副本未清除，单方生成影响 soundness | **代码层已修复，生产流程仍需落实** | RAII 擦除直接可达陷门副本；`keygen()` 停止在线生成 CRS；单方工具默认拒绝；加载时验证群、子群和幂链。生产仍必须使用经审计的多方 CRS，当前版本用于测试 |
| R9 | domain separator 和 `sid` 从 CRS seed 的 64-bit 前缀经 MT 派生 | **代码层已加固** | 完整 32-byte seed 经带用途标签的 SHAKE256 派生，不再截断到 64 bit，也不再使用 MT。若生产仪式格式允许，仍建议独立公开生成这些标签 |
| R10 | `HomConst.cpp` 存在未使用的弱随机引擎 | **已修复** | 删除无用随机引擎；常量初始化不再创建 RNG |
| R11 | 顶层构建未传播 TFHEpp 编译定义，存在 TU 配置不一致 | **已修复** | 顶层 CMake 传播所需 TFHEpp 编译定义和 OpenSSL 依赖 |
| R12 | 仓库残留包含旧 RNG 写法的备份源码 | **未进入构建，未作为运行路径** | 备份文件仍属于维护项；当前 CMake 目标不编译这些文件，正式发布打包应继续排除或删除 |

这里的“已修复”只表示当前第一方源码中的对应代码路径已经整改并通过回归。R8 的生产
多方可信设置是部署流程，不可能仅通过本地代码修改完成；R12 也没有被错误描述成已经
删除。

### 1.2 修复过程中另外发现的问题

| 问题 | 影响 | 处理结果 |
|---|---|---|
| TFHEpp BLAKE3 XOF refill 把元素数当成字节偏移，相邻缓冲区重叠 3584 字节 | 默认密钥生成、TLWE/TRLWE mask 和噪声底层随机流出现大面积确定性重复 | seek 改为按完整 4096-byte 缓冲区推进，并用临时独立探针对照连续 XOF 输出 |
| 64-bit Torus 噪声缩放表达式越界，实际噪声可能全部变为 0 | 压缩 LWE 公钥失去应有的误差项 | 使用 `std::ldexp` 在浮点域缩放，并增加 32/64-bit 回归 |
| ZK 可请求解压超过证明覆盖范围的槽位 | 返回未被证明约束的密文数据 | 同时限制 `P::n`、`params.k` 和实际证明范围 |
| C ABI 接受负长度、空指针，异常可能穿过 ABI | 越界读、崩溃或进程终止 | 统一校验参数、失败时清零输出长度，并在 ABI 内捕获异常 |
| 多个反序列化器信任声明长度或接受尾随数据 | 越界、超大分配、部分对象污染 | 增加上限、精确长度、剩余字节和整数溢出检查，使用临时对象提交 |
| 保存失败可能被误报为成功 | 产生不完整密钥集合但上层继续运行 | 传播内部返回值并检查打开、写入和 `flush` 状态 |
| 密钥生成/加载/保存线程中的异常可触发 `std::terminate` | 进程异常退出或对象部分初始化 | 线程捕获 `exception_ptr`，全部 `join` 后由主线程重抛 |
| 公钥文件中的字段可以覆盖全局 `CT_ZERO/CT_ONE` | 恶意公钥污染后续同态计算 | 旧字段只为格式兼容而消费，内容不再受信；常量在本地重建并加锁更新 |
| 密文包装对象缺少统一不变量 | 空向量、宽度不一致可能进入底层算子 | 统一检查非空、`W == data.size()`、宽度上限和操作数兼容性 |
| SHAKE256 的 EVP 返回值没有检查 | 加密提供者失败后可能继续使用无效输出 | 使用 RAII 管理上下文并检查全部 EVP 调用返回值 |

## 2. 修改后的随机源总体结构

当前代码把“需要公开且可确定性重建的数据”和“必须不可预测的秘密随机量”分成两条
独立路径：

```text
系统 CSPRNG
    │
    ├──取满 16 bytes──> PublicSeed
    │                       │
    │                       └──AES-128 counter stream──> 公开 a / mask
    │
    ├──SecureRng──> 公钥秘密误差 e
    ├──SecureRng──> 加密随机量 v、e0、e1、r_vec
    └──取满 32 bytes──> proof seed / CRS seed

b = a·s + e
压缩公钥 = (PublicSeed, b)
```

公开 seed 随公钥发布，任何人都可以用它重建 `a`，这是压缩公钥的正常设计。安全要求是
不能从该 seed 重建秘密误差 `e`；修改后的代码满足这一分离要求。

原生环境使用 OpenSSL `RAND_bytes` 取得系统安全随机字节，失败时抛出异常，不回退到弱
PRNG。Emscripten 环境循环调用由 `crypto.getRandomValues` 支持的 `random_device`，直到
填满请求长度，而不是只取一次 32-bit 输出。

## 3. 公开 seed 与秘密噪声分离

### 3.1 原来的压缩 LWE 公钥生成

原实现把公开 seed 写入公钥，然后用该 seed 为每个 `i` 创建 `mt19937_64`。公开矩阵行
`a_i` 和秘密噪声 `e_i` 来自同一个可重放随机流：

```cpp
std::uint64_t seed = std::random_device{}();
cpk.seed = seed;

for (size_t i = 0; i < M; ++i) {
    std::mt19937_64 rng(
        seed ^ (0x9E3779B97F4A7C15ULL + (i << 1))
    );

    for (size_t j = 0; j < P::n; ++j) {
        T aij = uniform(rng);
        inner += aij * s[j];
    }

    T e_i = gaussian(rng);
    cpk.b_list[i] = inner + e_i;
}
```

攻击者从公钥读取 seed 后可以执行完全相同的调用顺序，既重建 `a_i`，也精确重放
`e_i`，从而得到无噪声方程：

```text
b_i - e_i = <a_i, s>
```

这不是“seed 熵不够”这么简单，而是秘密噪声直接由公开数据确定，恢复噪声不需要暴力
搜索。

### 3.2 修改后的压缩 LWE 公钥生成

现在公开系数和秘密误差由两个互不派生的对象生成：

```cpp
SeededLwePublicKey generate_compressed_public_key(
    const SecretKey& sk,
    PublicSeed seed = secure_public_seed())
{
    SeededLwePublicKey cpk{};
    cpk.seed = seed;

    SecureRng noise_rng;
    PublicSeedStream public_stream(seed, 0x4c574541u); // "LWEA"

    for (size_t i = 0; i < M; ++i) {
        T inner = 0;
        for (size_t j = 0; j < P::n; ++j) {
            T aij = public_stream.next_uniform<T>();
            inner += aij * s[j];
        }

        T e_i = torus_from_normalized_noise<T>(gaussian(noise_rng));
        cpk.b_list[i] = inner + e_i;
    }
    return cpk;
}
```

两条路径为：

```text
PublicSeed ──AES-128 stream──> a_i

系统 CSPRNG ──SecureRng──────> e_i
```

公开 seed 仍然能够重建 `a_i`，但没有任何接口把它传给 `noise_rng`，因此不能再重放
`e_i`。

## 4. 低熵随机引擎替换

### 4.1 原问题的共同根因

原安全敏感入口多次使用：

```cpp
std::default_random_engine engine(std::random_device{}());
```

在原报告对应的 Emscripten/libc++ 环境中，`default_random_engine` 是状态空间约 31 bit
的线性同余引擎。一次 `random_device()` 只提供一个约 32-bit 值，之后产生的所有输出都
只是该短状态的确定性函数；输出再多也不会增加熵。

### 4.2 seeded RLWE 公钥噪声

原调用点：

```cpp
std::default_random_engine eng(std::random_device{}());
spk_lvl1 = generate_rlwe_pk_pair_seeded<lvl1param>(sk, eng);
spk_lvl2 = generate_rlwe_pk_pair_seeded<lvl2param>(sk, eng);
```

`generate_rlwe_pk_pair_seeded` 会用传入的 `rng` 生成公钥秘密误差多项式：

```cpp
spk.b = sample_gaussian_poly<P>(rng, P::α);
```

现在 lvl1 和 lvl2 调用点分别创建安全引擎：

```cpp
auto gen_spk_lvl1 = [&] {
    autohog::SecureRng eng;
    spk_lvl1 = generate_rlwe_pk_pair_seeded<lvl1param>(sk, eng);
};

auto gen_spk_lvl2 = [&] {
    autohog::SecureRng eng;
    spk_lvl2 = generate_rlwe_pk_pair_seeded<lvl2param>(sk, eng);
};
```

公钥生成函数继续保留外部 `URBG&` 参数，并实际消费调用点传入的 `SecureRng`：

```cpp
template<class P, class URBG>
SeededRlwePkPair<P> generate_rlwe_pk_pair_seeded(
    const SecretKey& sk,
    URBG& rng,
    PublicSeed seed = secure_public_seed())
{
    SeededRlwePkPair<P> spk{};
    spk.seed = seed;
    spk.b = sample_gaussian_poly<P>(rng, P::α); // 秘密 e

    for (size_t i = 0; i < P::k; ++i) {
        uniform_poly_from_seed<P>(a_i, spk.seed, i); // 公开 a_i
        spk.b += a_i * s_i;
    }
    return spk;
}
```

因此不是“保留了一个无用的兼容参数，然后在函数内部另建 RNG”。当前仓库内正式调用点
明确传入 `SecureRng`，函数内部直接使用它生成秘密噪声。

### 4.3 普通、V2 与 ZK 加密随机量

原普通、V2 和 ZK 入口分别创建 `default_random_engine`，该引擎直接生成：

```text
普通/V2：v、e0、e1
ZK：     r_vec、e0、e1
```

修改后三个入口均使用 `SecureRng`：

```cpp
trlwe_ptr encryptData(...) const {
    autohog::SecureRng eng;
    return encrypt_bits_trlwe_with_seeded_rlwe_pair(
        spk_lvl2, bits, eng);
}

trlwe_v2_ptr encryptDataV2(...) const {
    autohog::SecureRng eng;
    return encrypt_bits_trlwe_with_seeded_rlwe_pair(
        spk_lvl1, bits, eng);
}

Bytes encryptDataWithZk(...) const {
    autohog::SecureRng engine;
    return ZKPKE::encrypt_bits_with_zk(
        bits, spk_lvl2, engine, zk_params, metadata);
}
```

采样分布和 FHE 参数没有改变：仍然使用相同的高斯分布、`P::α`、二元随机量和明文编码；
只替换底层随机源。因此该修改不会改变既有噪声分析的参数，反而使实际采样更接近分析
中假设的独立随机样本。

### 4.4 Proof seed

原代码用多次 `random_device()` 逐字节填充 proof seed。虽然调用次数足够，但随机源逻辑
分散。现在统一为一次满长度安全请求：

```cpp
std::array<uint8_t, ZKProof::SEED_LEN> proof_seed{};
autohog::secure_random_bytes(
    proof_seed.data(), proof_seed.size());
```

## 5. 128-bit 公开 seed 与 AES 展开

### 5.1 公开 seed 的熵

原压缩 LWE seed：

```cpp
std::uint64_t seed = std::random_device{}();
```

原 seeded RLWE seed：

```cpp
std::uint64_t seed =
    std::mt19937_64(std::random_device{}())();
```

两者虽然存入 64-bit 字段，随机输入都只有一次约 32-bit `random_device()` 输出。现在公开
seed 固定为 16 字节：

```cpp
using PublicSeed = std::array<std::uint8_t, 16>;
static_assert(sizeof(PublicSeed) == 16);

inline PublicSeed secure_public_seed() {
    PublicSeed seed{};
    secure_random_bytes(seed.data(), seed.size());
    return seed;
}
```

这提供完整 128-bit seed 空间，满足本项目提出的 128-bit 熵要求。

### 5.2 原来的逐 `i` seed 派生

旧代码使用：

```cpp
seed ^ (0x9E3779B97F4A7C15ULL + (i << 1))
```

它的目的不是为每个 `i` 增加新熵，而是从一个公开根 seed 确定性得到不同的伪随机流，
从而只保存一个 seed 就能重建全部公开矩阵。不能为每个 `i` 单独调用
`std::random_device{}`，因为这些新 seed 必须全部保存或传输，否则解压端无法重建相同的
`a_i`；同时大量系统熵调用也没有必要。

真正的问题在于旧派生使用简单 XOR 和非密码学 `mt19937_64`，并且还把同一流继续用于
秘密噪声。

### 5.3 修改后的 `PublicSeedStream`

现在使用 128-bit seed 作为 AES-128 key，通过固定 counter block 生成连续伪随机流：

```text
counter_block = "AHG1" || BE32(domain) || BE64(block_index)

stream_block = AES-128-Encrypt(PublicSeed, counter_block)
```

用途域分开 LWE 与 RLWE：

```cpp
PublicSeedStream lwe_stream(seed, 0x4c574541u);  // "LWEA"
PublicSeedStream rlwe_stream(seed, 0x524c5745u); // "RLWE"
```

压缩 LWE 顺序消费同一条流：

```cpp
for (size_t i = 0; i < M; ++i)
    for (size_t j = 0; j < P::n; ++j)
        a[i][j] = stream.next_uniform<T>();
```

seeded RLWE 需要单独重建第 `i` 个多项式，因此定位到同一条流中互不重叠的区间：

```cpp
uint64_t byte_offset = i * P::n * sizeof(P::T);
PublicSeedStream stream(seed, 0x524c5745u, byte_offset);
```

每个 `i` 没有新建随机 seed，而是使用同一个 AES key 的不同 counter/字节区间。这个结构
参考 TFHE-rs 的 compression seed、AES-CTR mask generator 和独立 noise generator
设计，但 AutoHog 的 `"AHG1" || domain || block_index` 编码以及序列化格式是本项目自己
固定的协议，不宣称与 TFHE-rs 逐字节兼容。

### 5.4 `SecureRng` 实现

`SecureRng` 满足 `UniformRandomBitGenerator`，可以继续传给现有
`std::normal_distribution` 和 `std::bernoulli_distribution`：

```cpp
class SecureRng {
public:
    using result_type = std::uint64_t;

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() {
        return std::numeric_limits<result_type>::max();
    }

    result_type operator()() {
        if (buffer_exhausted()) {
            secure_random_bytes(buffer.data(), buffer.size());
        }
        return take_u64_from_buffer();
    }
};
```

当前实现缓冲 4096 字节 CSPRNG 输出，减少频繁系统调用；缓冲区在析构时显式擦除。它只
用于秘密随机量，不用于重建公开对象。

### 5.5 TFHEpp BLAKE3 随机流重叠修复

TFHEpp 默认启用 `USE_BLAKE3`，以 32-byte 秘密 seed 初始化 BLAKE3 keyed XOF，并把
XOF 输出缓存成 4096-byte 缓冲区。`BLAKE3PRNG<uint64_t>` 的一个缓冲区包含 512 个
`uint64_t`，原 refill 代码错误地把元素个数当成了 BLAKE3 `seek` 的字节偏移：

```cpp
static constexpr uint buffer_len = 4096;       // bytes
static constexpr uint buffer_result_len = 512; // uint64_t elements

// 修改前
blake3_hasher_finalize_seek(
    &hasher,
    counter * buffer_result_len,
    reinterpret_cast<uint8_t*>(buffer.data()),
    buffer_len);
```

BLAKE3 `seek` 的单位是字节，因此修改前实际读取的窗口为：

```text
buffer 0: XOF bytes [0,    4096)
buffer 1: XOF bytes [512,  4608)
overlap :           3584 bytes = 448 uint64_t values
```

独立探针在修复前确认第二批的 448 个重叠位置全部与第一批相等：

```text
repeated=448/448
```

由于项目默认构造 TFHEpp `SecretKey`，仅 `lvl0param::n = 711` 就会超过首次 512-output
refill 阈值；评估密钥、annihilate key 及直接 TFHEpp TLWE/TRLWE 加密还会继续消费该
随机流。因此这不是未调用的第三方死代码，也不是 seed 熵不足，而是高熵 seed 扩展后的
输出区间重叠。

修改后按完整字节缓冲区推进：

```cpp
// 修改后
blake3_hasher_finalize_seek(
    &hasher,
    counter * buffer_len,
    reinterpret_cast<uint8_t*>(buffer.data()),
    buffer_len);
```

现在依次读取 `[0, 4096)`、`[4096, 8192)` 等互不重叠区间。验证时使用临时独立探针和
固定 seed，先通过修复后的包装器取得 1024 个 `uint64_t`，再直接从 BLAKE3 XOF 的
`[0, 8192)` 取得参考字节，要求两个完整缓冲区逐字节一致。该验证不是依赖“随机碰撞
概率”的弱断言，而是确认包装器输出恰好等于连续 XOF 流；探针及编译产物在验证后删除，
没有作为额外测试文件加入项目。

## 6. 64-bit Torus 噪声缩放

原压缩 LWE 路径包含：

```cpp
constexpr int BITS = sizeof(T) * 8;
long long e_scaled = std::llround(
    gauss(rng) * static_cast<long double>(1ULL << BITS));
```

当 `T` 为 64 bit 时，`1ULL << 64` 的移位量等于类型宽度，属于未定义行为；后续尝试用
`exp2(64)` 后再收窄到 64-bit 整数，同样不能表示 `2^64`。在实际编译器上，该问题使
压缩 LWE 公钥噪声残差全部变成 0。

现在 LWE 和 RLWE 共用安全转换函数，在浮点域完成 `2^BITS` 缩放：

```cpp
template<class Torus>
Torus torus_from_normalized_noise(long double noise) {
    using Uns = std::make_unsigned_t<Torus>;
    constexpr int BITS = std::numeric_limits<Uns>::digits;

    const long double scaled = std::ldexp(noise, BITS);
    check_finite_and_long_long_range(scaled);

    const long long rounded = std::llround(scaled);
    return static_cast<Torus>(static_cast<Uns>(rounded));
}
```

真实 `lvl2param` 回归结果：

```text
修复前：samples=131264, nonzero_residuals=0
修复后：samples=131264, nonzero_residuals=131248
```

少量高斯样本量化为 0 是正常现象；关键安全条件是整个公钥噪声不能退化为全 0。

## 7. CRS 与 ZK 安全链

### 7.1 CRS seed 派生不再截断

原代码只取 32-byte `alpha_seed` 的前 8 字节，再用 `mt19937_64` 派生 `hash_*_ds` 和
`sid`：

```cpp
uint64_t prefix = 0;
memcpy(&prefix, alpha_seed.data(), sizeof(prefix));
std::mt19937_64 rng(prefix);
```

现在完整吸收全部 32 字节，并加入固定用途标签：

```text
input = "AUTOHOG-ZK-DOMAIN-SEPARATOR-V1" || alpha_seed[0..31]
derived = SHAKE256(input, required_bytes)
```

随后按固定顺序切分出 `hash_ds`、`hash_t_ds`、`hash_agg_ds`、`hash_lmap_ds`、
`hash_z_ds`、`hash_w_ds`、`hash_gamma_ds` 和 `sid`。因此派生强度不再受 64-bit 前缀限制，
也不再依赖 Mersenne Twister。

### 7.2 陷门副本擦除

原实现只擦除了输入 seed，但 `alpha`、`alpha_fr` 及 OpenMP 工作线程中的
`start_power/start_scalar` 会保留在内存或线程栈上。现在使用不可复制的 RAII 包装器：

```cpp
template<class T>
class ScrubbedValue {
public:
    ~ScrubbedValue() {
        autohog::secure_erase(&value_, sizeof(value_));
    }
    T& get() { return value_; }
private:
    T value_{};
};
```

顶层和每个并行 chunk 都使用该包装器，因此正常返回和异常展开都会擦除直接可达副本。

### 7.3 禁止把在线单方 CRS 当作生产 CRS

`FheKey::keygen()` 现在只生成 FHE 密钥，不再同时生成单方 CRS。仓库 `gen_crs` 工具默认
拒绝单方生成，只有显式指定危险测试开关时才允许生成本地回归制品。

生产加载必须提供独立可信渠道获得的 SHA-256：

```cpp
key.load(
    eval_key_file,
    annihilate_key_file,
    public_key_file,
    public_key_v2_file,
    crs_file,
    trusted_crs_sha256);
```

无摘要的高层文件、内存和旧 C ABI 入口在请求加载 CRS 时失败关闭。文件认证路径只打开
CRS 一次，对同一内存缓冲区计算摘要并解析，避免“先哈希、再重新打开”产生文件替换窗口。
当前测试中仍使用单方生成的CRS。

### 7.4 CRS 结构验证

原验证只检查群元素是否合法，不能发现合法点被换序或不属于同一幂链。现在加载过程依次
检查：

1. 参数尺寸和受检算术；
2. G1/G2 点编码、群和子群；
3. `g_list[n]` 的规定洞位；
4. G1/G2 起点一致性；
5. G1 和 G2 相邻幂关系；
6. 洞位两侧的桥接关系。

相邻关系使用从完整 CRS transcript 经 SHAKE256 派生的 128-bit 批处理系数完成配对检查。
只有摘要和全部结构验证均通过后，才设置 `integrity_checked = true`。

代码擦除和幂链验证不能证明生成者已经销毁陷门。生产必须使用与当前
BLST/BLS12-381 参数及本项目索引布局兼容的多方仪式 CRS，并通过独立渠道固定其摘要。
原 PDF 中针对另一曲线提出的 BN254 `.ptau` 不能直接导入当前 BLS12-381 实现。

## 8. 外部输入、保存加载与运行时边界

这些问题不是原报告 R1-R9 的 RNG 根因，但与修复后的密钥和证明路径直接相关，已经一并
处理。

### 8.1 ZK 输出范围

原验证只限制 `num_bits <= P::n`，没有要求输出位数位于证明覆盖的 `params.k` 内。现在：

```cpp
if (num_bits == 0 ||
    num_bits > P::n ||
    num_bits > params.k) {
    return false;
}
```

加密端也拒绝超过 CRS 证明容量的输入。

### 8.2 C ABI 与批量格式

所有导出入口现在先检查 key、输入指针、输出指针和非负长度，失败时把 `*out_len` 置零；
整个 C ABI 捕获 C++ 异常并转换为失败返回，不允许异常穿过语言边界。

批量整数格式不再通过未对齐指针直接读取，也不再忽略 `in_len`。解析器维护明确的
`remaining`，使用大端逐字节读取，并限制 metadata、元素数量和输出长度，拒绝截断与
尾随数据。

### 8.3 反序列化与固定布局文件

普通 TLWE 字节输入必须非空且为 `sizeof(tlwe_t)` 的精确整数倍：

```cpp
if (bytes.empty() ||
    bytes.size() % sizeof(tlwe_t) != 0) {
    throw std::invalid_argument("invalid TLWE ciphertext byte length");
}
```

压缩密文、公钥、TLWE、TRLWE、RLWE 和 seeded RLWE 公钥均检查精确布局长度。proof 和
metadata 各有 1 MiB 上限。所有加载先解析到临时对象，完全成功后才更新调用者对象，
避免截断文件留下部分状态。

### 8.4 保存错误传播与并发异常

原部分保存入口忽略内部返回值并直接返回成功：

```cpp
saveSecretKey(sk, ofs);
ofs.flush();
return true;
```

现在检查内部序列化、流状态和 `flush`：

```cpp
if (!saveSecretKey(sk, ofs)) return false;
ofs.flush();
if (!ofs) return false;
return true;
```

密钥生成、加载和保存线程在内部捕获 `std::exception_ptr`，主线程等待全部线程结束后重新
抛出第一个错误。任一必需文件失败都不会被当作成功，也不会让异常逃出线程函数触发
`std::terminate`。

### 8.5 全局同态常量

旧 lvl2 seeded 公钥格式尾部保存 `CT_ZERO/CT_ONE`，加载代码曾直接信任文件内容并覆盖
进程全局变量。现在仍消费这些字段以维持旧布局长度，但不使用其值，而是重建确定性
trivial 常量：

```cpp
TLWE<lvl2param> zero{};
zero[lvl2param::n] = -lvl2param::μ;

TLWE<lvl2param> one{};
one[lvl2param::n] = lvl2param::μ;
```

常量从静态初始化开始即有效，更新集中到 `init_constants` 并使用互斥锁。这样切断了恶意
公钥污染全局状态的路径；不过它们仍是进程级对象，如果系统需要同一进程并行服务多个
独立密钥，后续仍应迁移为每个 `FheKey` 实例自己的状态。

### 8.6 密文包装对象与 SHAKE 错误

公开算术、比较、位运算、压缩、解压和解密入口统一检查：

```text
对象指针非空
data 非空
W == data.size()
1 <= W <= lvl2param::n
二元操作的左右宽度相等
select 条件宽度为 1
```

同时修复空批次访问首元素、反向循环的有符号下溢和整数向量反序列化长度不精确问题。

SHAKE256 辅助函数现在使用 RAII 管理 `EVP_MD_CTX`，并要求
`EVP_DigestInit_ex`、`EVP_DigestUpdate`、`EVP_DigestFinalXOF` 全部返回 1；任何失败都会
立即抛错，不继续使用无效挑战数据。

## 9. 文件格式和部署影响

### 9.1 公钥格式

公开 seed 字段从 8 字节变为 16 字节：

```text
旧格式：[8-byte seed][public-key body]
新格式：[16-byte seed][public-key body]
```

同时，公开 `a` 的确定性解释从 MT/XOR 派生改为 AES-128 计数流。即使某个中间版本已经
使用 16-byte seed，只要公开展开算法不同，也不能和当前版本混用。

升级必须重新生成 compressed LWE 公钥和 lvl1/lvl2 seeded RLWE 公钥，并同步替换客户端、
服务端和测试夹具。若历史公钥由公开 seed 同时派生秘密噪声的代码生成，应按对应私钥
可能泄露处理，轮换私钥而不是只替换公钥文件。

### 9.2 CRS

CRS 的 domain separator 派生算法、加载信任锚和结构验证已经改变。证明端与验证端必须
使用同一份通过新验证器的 CRS，并从独立可信渠道固定 SHA-256。测试用单方 CRS 不能
替代生产多方仪式制品。

### 9.3 测试数据

更换安全随机引擎会改变随机密文和 proof 字节，因此依赖固定随机输出的黄金值需要重新
建立；确定性公开 seed 测试则使用固定向量验证 AES 流的一致性。

## 10. 修改文件汇总

| 文件 | 主要修改 |
|---|---|
| `include/SecureRandom.h` | 安全字节入口、128-bit `PublicSeed`、AES 公开流、缓冲式 `SecureRng`、安全擦除 |
| `include/HomPKE.h` | 压缩 LWE/RLWE 的公开流与秘密噪声分离、AES 展开、Torus 噪声转换 |
| `include/ZKPKE.h` | ZK 加密随机量、proof seed 和证明范围检查 |
| `include/fhe_key.h` | 安全 RNG 调用点、CRS 认证入口、生成/加载/保存异常传播、密文对象不变量 |
| `include/SaveFile.h` | 密钥、密文、proof、metadata 的有界保存加载与精确文件布局 |
| `include/ZKProof.h` | CRS 认证加载与完整性状态接口 |
| `include/fhe_key_api.h` | 带固定 CRS SHA-256 的认证 C ABI |
| `src/fhe_key.cpp` | 公钥/密文序列化边界、文件长度、全局常量污染修复 |
| `src/fhe_key_api.cpp` | C ABI 指针长度检查、批量解析、异常边界和输出长度复位 |
| `src/ZKProof.cpp` | SHAKE 检查、CRS seed 完整派生、陷门擦除、摘要绑定和幂链验证 |
| `src/HomConst.cpp` | 删除无用 RNG，确定性常量初始化和加锁更新 |
| `CMakeLists.txt` | TFHEpp 编译定义和 OpenSSL 依赖传播 |
| `thirdparties/TFHEpp/include/BLAKE3PRNG.hpp` | 本地上游补丁：BLAKE3 refill seek 按 4096-byte 缓冲区推进，消除 3584-byte 重叠 |
| `Coprocessor-operator/gen_crs.cpp` | 单方生成默认拒绝、测试危险开关、摘要输出和验证重载 |
| `Coprocessor-operator/test_public_seed_stream.cpp` | AES 公开流固定向量、连续区间和域分离测试 |
| `Coprocessor-operator/test_PKE.cpp` | 标准/压缩 LWE、RLWE 和 32/64-bit Torus 缩放回归 |
| `Coprocessor-operator/test_with_zk.cpp` | ZK 全流程、篡改拒绝和证明范围测试 |
| `Coprocessor-operator/test_security_boundaries.cpp` | C ABI、反序列化、畸形包装器、伪造常量和保存失败负向测试 |


## 11. 构建与测试结果


| 项目 | 结果 |
|---|---|
| `homops`、`fhe-api` | BUILD PASS |
| `test_public_seed_stream` | PASS |
| `test_PKE` | PASS；标准/压缩 LWE、标准/seeded RLWE 全部通过 |
| 32/64-bit Torus 噪声缩放 | PASS |
| `test_security_boundaries` | PASS |
| `test_fhe_key_api` | PASS |
| `test_pack_unpack` | PASS |
| `test_with_zk 32` | PASS；加密、证明、验证、序列化、解密和 bootstrap 通过 |
| 请求证明范围之外的输出 | 正确拒绝 |
| CRS 幂链、proof、metadata、公钥篡改 | 正确拒绝 |
| 截断/尾随文件、空向量、异宽操作、保存失败 | 正确拒绝 |

测试使用的临时单方 CRS、测试密钥、密文和构建目录已在验证结束后清理。项目当前没有把
这些可执行测试完整注册进 CTest，因此上述结果来自直接构建并运行目标文件。

## 12. 当前仍需保留的安全边界

1. 生产必须使用经独立验证的多方 BLS12-381 CRS，并通过可信渠道固定最终摘要。代码中的
   单方 CRS 测试开关只能用于回归。
2. `generate_rlwe_pk_pair_seeded` 保留通用 `URBG&` 模板。仓库内正式调用点全部传入
   `SecureRng`；仓库外调用者如果直接实例化该模板，也必须满足同一要求。
3. `CT_ZERO/CT_ONE` 已不再接受公钥文件内容，但仍为进程全局对象。多密钥、多租户并发
   场景应进一步改为实例状态。
4. 秘密密钥文件仍需要由部署方配置严格权限、原子替换、备份策略和安全内存处理。这是
   密钥存储加固，不是本次 RNG 熵源修复本身。
5. 含旧实现的备份源码不参与当前构建，但发布打包和后续维护应明确排除或删除，防止误用。
6. TFHEpp 时必须保留该补丁，直到上游版本包含等价修复，并重新核对跨 refill 输出。

## 13. 最终结论

原 `RNG_AUDIT.pdf` 指出的核心密码学随机数问题已经在当前第一方运行路径中得到修复：

- 公开 seed 不再派生秘密噪声；
- 公钥、普通加密、V2 加密和 ZK 加密的秘密随机量不再来自 31-bit LCG；
- 公开 seed 具有完整 128-bit 熵，并使用密码学 AES 计数流展开；
- 短状态随机流重叠的原根因已经消除；
- CRS seed 不再截断，直接可达陷门副本得到擦除，加载端具备摘要和幂链验证；
- 64-bit Torus 公钥噪声不再退化为 0；
- 与密钥、密文和证明相关的外部输入、保存加载和异常边界已经补齐。

需要继续由系统设计和部署保证的事项是生产多方 CRS、全局同态常量的实例隔离、秘密
密钥存储策略以及备份源码治理。它们在本报告中被明确列为边界，没有被计入已经完成的
代码修复。
