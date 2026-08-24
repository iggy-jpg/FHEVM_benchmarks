# 密文并发一致性修复报告：固定 FFTW Wisdom

## 1. 结论

本次修复不再依赖“清除密文系数最低 4 bit”来掩盖并发执行产生的差异，而是固定 FFTW 的规划结果（plan），让所有工作线程复用由同一份 Wisdom 约束的 FFT 实现路径。

最终结果如下：

- 在关闭明文折叠和关闭低 4 bit 归零后，并行度 2、4、8、16 的原始密文均完全一致。
- 密文乘法和密文除法已包含在严格一致性测试中，所有并行度下均为 0 个 raw mismatch。
- 所有测试结果的明文解密值均正确。
- 汇总 22 个主要密文算子的平均耗时后，新方案相对原始 FFTW `MEASURE` 的变化为 **+0.13%**，基本保持原有性能。
- FFTW plan 只在线程首次使用 FFT 时创建，随后被该线程反复复用；不会在每次调用算子前重新创建。

本次最终严格一致性结果位于：

```text
results/consensus_raw_20260805_211747
```

---

## 2. 原有的低 4 bit 归零机制

原有代码提供了 `sanitize_tlwe` 处理，对 TLWE 密文的每个系数执行：

```cpp
val &= ~0xFULL;
```

它会将每个系数最低 4 bit 清零，再比较处理后的密文。

这种机制的出发点是：如果不同线程只在最低几位出现微小数值差异，那么丢弃这些位后，密文可以被判定为一致。由于这些低位通常位于噪声容限内，解密结果可能仍然相同。

但它存在三个问题：

1. 它只处理比较结果，没有消除差异产生的根源。
2. 它会丢弃密文系数的低位信息，严格意义上不能证明原始密文字节一致。
3. 乘法、除法等复杂算子的误差会经过多级同态计算传播，不能仅靠固定清除 4 bit 保证所有情况都安全。

因此，本次验收使用：

```text
--disable-low4-sanitize
```

测试直接比较未经归零处理的原始密文。

---

## 3. 不一致的原因

FFTW 的 `MEASURE` 规划会在创建 plan 时实际试跑多个候选算法，并选择当时测得较快的实现。并发线程分别创建 plan 时，机器负载、缓存状态和测量抖动可能导致线程选中不同的 FFT 算法或执行路径。

这些实现路径在数学意义上等价，但浮点运算的结合顺序、舍入误差和 SIMD 执行顺序可能不同。FFT 结果再转换回 Torus 整数时，微小浮点差异可能落到密文系数的最低若干位，最终表现为：

- 解密明文一致；
- 原始密文的部分系数不完全一致；
- 清除最低 4 bit 后差异消失。

因此，问题的核心不是同态算子逻辑随机变化，而是不同线程使用了不同的 FFTW plan。

---

## 4. 新方案：一次校准，固定 Wisdom，线程内复用 plan

新方案保留 `MEASURE` 对性能较好的算法搜索，但将搜索结果固化为一份 FFTW Wisdom，并要求后续线程只能使用 Wisdom 中记录的 plan。

### 4.1 执行流程

1. 进程读取 `AUTOHOG_FFTW_WISDOM_FILE` 指定的 Wisdom 文件。
2. 如果文件不存在或与当前环境不兼容，第一个进入 FFTW 初始化的线程执行一次 `MEASURE` 校准。
3. 校准完成后，将 plan 信息导出到 Wisdom 文件。
4. 后续线程导入同一份 Wisdom，并以 `FFTW_WISDOM_ONLY` 创建 plan。
5. 每个线程的 FFTW processor 和 plan 使用 `thread_local` 保存，在该线程后续的所有算子中复用。

### 4.2 并发安全

FFTW 的 plan 创建、销毁以及 Wisdom 导入、导出由同一个全局互斥锁保护，避免多个线程同时修改 FFTW planner 的全局状态。

FFT 实际执行不持有该锁，各个线程使用自己的 plan 并行执行，因此不会把算子执行串行化。

### 4.3 与低 4 bit 归零的本质区别

| 方案 | 处理位置 | 是否保留原始密文 | 是否消除根因 | 对算子执行性能的影响 |
|---|---|---:|---:|---:|
| 低 4 bit 归零 | FFT 计算完成后 | 否 | 否 | 很小，但只是容忍差异 |
| 固定 Wisdom | FFTW plan 创建阶段 | 是 | 是 | 主要密文算子平均 +0.13% |

新方案不是把差异截掉，而是确保所有线程选择相同的 FFT 算法路径，使差异不再产生。

---

## 5. Wisdom 何时创建

构建项目时不会创建 Wisdom。它在程序运行后、线程第一次真正使用 FFT 时创建或加载。

在 `test_concurrency.sh` 中，流程为：

1. 脚本创建临时工作目录和临时 Wisdom 路径。
2. 先运行外部并行度 1。
3. 第一次 FFT 初始化时生成 Wisdom；同一进程的其他线程读取并复用该结果。
4. 外部并行度 2、4、8、16 的测试继续复用同一份 Wisdom。
5. 测试结束后，Wisdom 被复制到最终结果目录，临时工作目录被清理。

所以，每次完整执行测试脚本只需要做一次校准，不是每个算子一次，也不是每次 FFT 运算一次。

在长期运行的服务中，可以把 `AUTOHOG_FFTW_WISDOM_FILE` 指向持久化文件。只要 CPU、FFTW 版本和关键构建选项未变化，后续进程就可以直接导入；环境变化后应重新校准。

---

## 6. 测试代码和数据调整

### 6.1 一致性测试

`Coprocessor-operator/test_consensus.cpp` 当前支持：

- 分别统计 raw、低 4 bit 归零后以及明文结果的不一致数量；
- 固定常量密文缓存，避免常量生成引入额外随机性；
- `mul_pal` 和 `div_pal` 真正使用密文输入，不再只覆盖明文折叠路径；
- 固定总线程预算为 16：外部并行度 1、2、4、8、16 分别对应每个 worker 使用 16、8、4、2、1 个算子线程；
- 每个 worker 使用独立临时目录，结束后自动清理。

`test_concurrency.sh` 负责顺序运行各个并行度、复用相同 Wisdom、汇总结果并检查失败状态。

### 6.2 严格一致性工作负载

`data/consensus_operations_consistency.txt` 为每个受支持算子提供固定输入，覆盖：

```text
add、sub、mul、div、and、or、not、ge、lt、select、常量和上传
```

测试时同时关闭明文折叠与低 4 bit 归零，确保乘法和除法也进入真实密文计算路径，并按原始密文比较。

### 6.3 性能测试

`Coprocessor-operator/test_PERF_ALL.cpp` 使用固定随机种子，并在正式计时前执行一次不计入结果的预热，减少随机输入和首次调用对结果的影响。

---

## 7. 密文一致性实测结果

测试条件：

- `plain_fold=0`
- `low4_sanitize=0`
- `shared_key=1`
- 总线程预算为 16
- 所有并行度复用同一份 Wisdom
- Wisdom SHA-256：`0b3c2bf5416c531a1eb484397b9c1678454f242c0ccc78c46d7d9c71ba051cb6`

| 外部并行度 | 比较次数 | raw mismatch | sanitized mismatch | plaintext mismatch |
|---:|---:|---:|---:|---:|
| 2 | 32 | 0 | 0 | 0 |
| 4 | 64 | 0 | 0 | 0 |
| 8 | 128 | 0 | 0 | 0 |
| 16 | 256 | 0 | 0 | 0 |

其中密文乘法和密文除法在 2、4、8、16 四种并行度下的 raw mismatch 均为 0。所有算子的预期明文检查也均为 PASS。

这说明一致性是在未使用低 4 bit 归零的条件下取得的。

---

## 8. 算子速度对比

### 8.1 测试方法

- FFTW：3.3.10-5；
- OpenMP 线程数：8；
- 算子宽度：8；
- 分别构建原始 `MEASURE` 和新 Wisdom 两套二进制；
- Wisdom 的首次校准不计入算子耗时；
- 每种方案启动 3 个独立进程，每个进程预热后记录 3 个 epoch；
- 每个算子、每种方案共 9 个有效样本；
- 全部 34 个算子的正确性检查通过。

### 8.2 主要密文算子

单位为毫秒，正数表示新方案更慢，负数表示新方案更快。

| 算子 | 原始 MEASURE | 新 Wisdom | 变化 |
|---|---:|---:|---:|
| add_s | 91.740 | 95.300 | +3.88% |
| add_u | 97.463 | 94.883 | -2.65% |
| and | 24.193 | 24.887 | +2.87% |
| div | 950.943 | 954.857 | +0.41% |
| div_newton | 4859.600 | 4844.437 | -0.31% |
| eq | 46.607 | 47.153 | +1.17% |
| ge_s | 47.360 | 47.410 | +0.11% |
| ge_u | 46.917 | 47.603 | +1.46% |
| gt_s | 47.397 | 48.240 | +1.78% |
| gt_u | 47.150 | 48.653 | +3.19% |
| le_s | 47.303 | 46.853 | -0.95% |
| le_u | 47.450 | 49.313 | +3.93% |
| lt_s | 46.983 | 47.137 | +0.33% |
| lt_u | 46.943 | 47.063 | +0.26% |
| mul | 323.753 | 329.310 | +1.72% |
| mul_low | 215.280 | 218.690 | +1.58% |
| ne | 47.040 | 47.497 | +0.97% |
| or | 24.130 | 24.867 | +3.05% |
| select | 32.867 | 34.637 | +5.39% |
| sub_s | 92.157 | 93.660 | +1.63% |
| sub_u | 95.280 | 94.990 | -0.30% |
| xor | 24.170 | 25.143 | +4.03% |

22 个主要密文算子各自平均耗时之和：

| 方案 | 汇总耗时 |
|---|---:|
| 原始 MEASURE | 7302.726 ms |
| 新 Wisdom | 7312.583 ms |
| 变化 | **+0.13%** |

### 8.3 密文与明文混合算子

| 算子 | 原始 MEASURE | 新 Wisdom | 变化 |
|---|---:|---:|---:|
| div_plain | 392.720 | 394.380 | +0.42% |
| div_plain_q | 181.153 | 181.773 | +0.34% |
| eq_plain | 22.247 | 21.853 | -1.77% |
| gt_plain | 45.027 | 45.833 | +1.79% |
| mul_plain | 181.840 | 186.470 | +2.55% |
| mul_plain_low | 113.647 | 114.477 | +0.73% |

6 个混合算子各自平均耗时之和从 936.634 ms 变为 944.786 ms，变化为 **+0.87%**。

### 8.4 极短算子说明

`and_plain`、`cast`、`not`、`or_plain`、`select_plain` 和 `xor_plain` 的单次耗时约为 0.01～0.06 ms。其绝对差异只有几十微秒，百分比容易被计时粒度、调度和缓存噪声放大，不适合用百分比判断 FFTW 策略的性能。

总体上，新方案没有改变 FFT 执行阶段的算法，只固定了 plan 的选择；因此校准完成后，主要算子性能与原始 `MEASURE` 基本一致。

---

## 9. 完整运行指令

以下命令均在 `AutoHog-Coprocessor` 项目根目录执行。

### 9.1 安装依赖和拉取子模块

```bash
apt update
apt install -y \
  build-essential cmake ninja-build \
  libfftw3-dev libssl-dev libboost-dev python3

git submodule update --init --recursive
```

### 9.2 构建一致性和性能测试程序

```bash
cmake -S . -B build-consensus -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_FFTW3=ON \
  -DTFHEPP_FFTW_PLAN_MODE=WISDOM \
  -DTFHEPP_ENABLE_IPO=OFF

cmake --build build-consensus \
  --target test_consensus test_PERF_ALL \
  -j "$(nproc)"
```

### 9.3 运行严格密文一致性测试

```bash
bash test_concurrency.sh \
  --computations data/consensus_operations_consistency.txt \
  --disable-plain-fold \
  --disable-low4-sanitize \
  --shared-key
```

脚本默认测试外部并行度 1、2、4、8、16。并行度 1 作为基准，2、4、8、16 与基准密文比较。

如果密钥不存在，`test_consensus` 会自动生成。需要主动重新生成全部测试密钥时使用：

```bash
bash test_concurrency.sh \
  --computations data/consensus_operations_consistency.txt \
  --disable-plain-fold \
  --disable-low4-sanitize \
  --shared-key \
  --regenerate-keys
```

### 9.4 运行原始完整工作负载

按原有行为运行 `data/computations.txt`：

```bash
bash test_concurrency.sh \
  --computations data/computations.txt \
  --shared-key
```

对完整工作负载执行最严格的原始密文检查：

```bash
bash test_concurrency.sh \
  --computations data/computations.txt \
  --disable-plain-fold \
  --disable-low4-sanitize \
  --shared-key
```

`data/computations.txt` 是较大的完整计算图；关闭明文折叠后会让更多节点进入真实密文计算，运行时间可能从数小时到更久。快速回归应优先使用 `data/consensus_operations_consistency.txt`。

---

## 10. 修复工作总结

本次为修复并验证密文一致性完成了以下工作：

1. 分离 raw、低 4 bit 归零后和明文三个层次的一致性检查。
2. 增加真实密文乘法、密文除法测试数据，避免被明文折叠绕过。
3. 确认低 4 bit 归零只能容忍低位差异，不能解决 FFTW plan 不一致的根因。
4. 在底层 FFTW processor 中加入 Wisdom 导入、校准、导出和严格复用逻辑。
5. 对 FFTW planner 的全局操作加锁，同时保留 FFT 执行阶段的线程并行能力。
6. 使用线程局部 processor，让 plan 每线程只初始化一次并持续复用。
7. 调整并发测试脚本，统一 Wisdom、线程预算、临时目录和结果归档。
8. 在关闭低 4 bit 归零后完成 2、4、8、16 并行度的严格原始密文一致性测试。
9. 使用独立构建和多进程、多 epoch 方法，对原始 `MEASURE` 与新 Wisdom 方案进行算子性能对比。

最终方案在不修改密文低位的前提下实现了并发密文一致性，并将主要密文算子的汇总性能变化控制在 **+0.13%**。
