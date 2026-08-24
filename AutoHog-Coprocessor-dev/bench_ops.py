"""
✅ bench_ops.py (多次测试取平均版)
----------------------------------------
用于批量测试 Coprocessor-operator 的算子性能。

支持算子：
  ADD, SUB, AND, EQ, GT, NOT, SELECT, OR, XOR, LT, GE, LE

支持模式：
  - 每个算子会被调用其对应的可执行文件（test_XXX）
  - 支持多次运行取平均、计算标准差
  - 脚本会自动解析输出中的 "AVG: <number> ms" 行，提取平均耗时

新增功能：
  ✅ --repeat 参数：每个测试重复N次（默认10次）
  ✅ 统计信息：平均值、标准差、最小值、最大值
  ✅ --warmup 参数：预热运行（不计入统计）

使用示例：
  python3 bench_ops.py \
    --ops SELECT \
    --bits 64 \
    --repeat 20 \
    --warmup 3 \
    --csv results/select_stable.csv
"""

import argparse
import csv
import re
import subprocess
import sys
import statistics
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

# -----------------------------

# -----------------------------
DEFAULT_BINARIES = {
    "ADD":    "./build/Coprocessor-operator/test_ADD",
    "SUB":    "./build/Coprocessor-operator/test_SUB",
    "AND":    "./build/Coprocessor-operator/test_AND",
    "EQ":     "./build/Coprocessor-operator/test_EQ",
    "GT":     "./build/Coprocessor-operator/test_GT",
    "NOT":    "./build/Coprocessor-operator/test_NOT",
    "SELECT": "./build/Coprocessor-operator/test_SELECT",
    "OR":     "./build/Coprocessor-operator/test_OR",
    "XOR":    "./build/Coprocessor-operator/test_XOR",
    "LT":     "./build/Coprocessor-operator/test_LT",
    "GE":     "./build/Coprocessor-operator/test_GE",
    "LE":     "./build/Coprocessor-operator/test_LE",
}

# -----------------------------

# -----------------------------
REGEX_PATTERNS = {
    "ADD": {
        "SER": re.compile(r"HomAdd\s*\(serial\)\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomAdd\s*\(parallel\)\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "SUB": {
        "SER": re.compile(r"SER\s*avg\s*time\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"PAR\s*avg\s*time\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "AND": {
        "AVG": re.compile(r"HomAND\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "EQ": {
        "SER": re.compile(r"HomEQ\s*\[SER\][^\n]*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomEQ\s*\[PAR\][^\n]*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "GT": {
        "SER": re.compile(r"HomGT\s*\[SER\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomGT\s*\[PAR\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "NOT": {
        "AVG": re.compile(r"Average\s*time\s*over\s*\d+\s*epoch\(s\)\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "SELECT": {
        "AVG": re.compile(r"HomSelect\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "OR": {
        "AVG": re.compile(r"HomOR\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "XOR": {
        "AVG": re.compile(r"Average\s*time\s*over\s*\d+\s*epoch\(s\)\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "LT": {
        "SER": re.compile(r"HomLT\s*\[SER\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomLT\s*\[PAR\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "GE": {
        "SER": re.compile(r"HomGE\s*\[SER\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomGE\s*\[PAR\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
    "LE": {
        "SER": re.compile(r"HomLE\s*\[SER\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
        "PAR": re.compile(r"HomLE\s*\[PAR\]\s*execution\s*AVG\s*:\s*([0-9]*\.?[0-9]+)\s*ms", re.IGNORECASE),
    },
}

@dataclass
class BenchResult:
    op: str
    bits: str
    epochs: Optional[int]
    metric: str
    mean_ms: float
    stdev_ms: float
    min_ms: float
    max_ms: float
    repeat: int
    cmd: str
    rc: int
    notes: str = ""
    raw_values: List[float] = field(default_factory=list)

def _build_cmd(op: str, binary: str, bits: str, epochs: Optional[int], detail_flag: Optional[str]) -> List[str]:
    """
    根据各自用法拼装命令行参数。
    """
    args = [binary]
    args.append(bits)
    if epochs is not None:
        args.append(str(epochs))
    if detail_flag in ("--no-epoch-detail", "--show-epochs"):
        args.append(detail_flag)
    return args

def _run_once(cmd: List[str]) -> Tuple[int, str, str]:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return proc.returncode, proc.stdout, proc.stderr

def _parse_output(op: str, stdout: str) -> Dict[str, float]:
    patterns = REGEX_PATTERNS.get(op.upper())
    if patterns is None:
        raise ValueError(f"Unknown op for parsing: {op}")
    results: Dict[str, float] = {}
    for metric, rx in patterns.items():
        m = rx.search(stdout)
        if m:
            try:
                results[metric] = float(m.group(1))
            except ValueError:
                pass
    return results

def run_test(
    op: str,
    bits: Union[int, str],
    epochs: Optional[int] = None,
    binaries: Optional[Dict[str, str]] = None,
    detail_flag: Optional[str] = None,
    repeat: int = 10,
    warmup: int = 2,
    verbose: bool = False,
) -> List[BenchResult]:
    """
    运行单个测试，支持多次重复和预热。
    
    Args:
        op: 操作名称
        bits: 位宽
        epochs: 测试轮数
        binaries: 自定义二进制路径
        detail_flag: 详细输出标志
        repeat: 重复次数
        warmup: 预热次数
        verbose: 是否显示每次运行结果
    """
    op = op.upper()
    bins = dict(DEFAULT_BINARIES)
    if binaries:
        bins.update({k.upper(): v for k, v in binaries.items()})
    if op not in bins:
        raise ValueError(f"Unsupported op: {op}")

    if isinstance(bits, int):
        bits_str = str(bits)
    elif isinstance(bits, list):
        bits_str = ",".join(map(str, bits))
    else:
        bits_str = str(bits)

    cmd = _build_cmd(op, bins[op], bits_str, epochs, detail_flag)
    

    if warmup > 0 and verbose:
        print(f"🔥 Warming up {warmup} times...", file=sys.stderr)
    for i in range(warmup):
        rc, out, err = _run_once(cmd)
        if verbose:
            print(f"  Warmup {i+1}/{warmup}: done", file=sys.stderr)
    

    all_metrics: Dict[str, List[float]] = {}
    failed_runs = 0
    
    if verbose:
        print(f"🏃 Running {repeat} times...", file=sys.stderr)
    
    for i in range(repeat):
        rc, out, err = _run_once(cmd)
        
        if rc != 0:
            failed_runs += 1
            if verbose:
                print(f"  Run {i+1}/{repeat}: FAILED (rc={rc})", file=sys.stderr)
            continue
            
        metrics = _parse_output(op, out)
        
        if not metrics:
            failed_runs += 1
            if verbose:
                print(f"  Run {i+1}/{repeat}: No metrics parsed", file=sys.stderr)
            continue
        

        for metric_name, value in metrics.items():
            if metric_name not in all_metrics:
                all_metrics[metric_name] = []
            all_metrics[metric_name].append(value)
        
        if verbose:
            print(f"  Run {i+1}/{repeat}: {metrics}", file=sys.stderr)
    

    results: List[BenchResult] = []
    
    if not all_metrics:
        sys.stderr.write(f"[WARN] 未解析到任何有效结果，op={op} bits={bits_str}\n")
        results.append(
            BenchResult(
                op=op,
                bits=bits_str,
                epochs=epochs,
                metric="N/A",
                mean_ms=float("nan"),
                stdev_ms=float("nan"),
                min_ms=float("nan"),
                max_ms=float("nan"),
                repeat=repeat - failed_runs,
                cmd=" ".join(cmd),
                rc=rc,
                notes=f"Failed runs: {failed_runs}/{repeat}",
                raw_values=[],
            )
        )
    else:
        for metric_name, values in all_metrics.items():
            if len(values) < 2:
                stdev = 0.0
            else:
                stdev = statistics.stdev(values)
            
            results.append(
                BenchResult(
                    op=op,
                    bits=bits_str,
                    epochs=epochs,
                    metric=metric_name,
                    mean_ms=statistics.mean(values),
                    stdev_ms=stdev,
                    min_ms=min(values),
                    max_ms=max(values),
                    repeat=len(values),
                    cmd=" ".join(cmd),
                    rc=0,
                    notes=f"Success: {len(values)}/{repeat}, Failed: {failed_runs}/{repeat}",
                    raw_values=values,
                )
            )
    
    return results

def run_benchmarks(
    ops: List[str],
    bits: Union[int, str, List[int]],
    epochs: Optional[int],
    binaries: Optional[Dict[str, str]] = None,
    detail_flag: Optional[str] = None,
    repeat: int = 10,
    warmup: int = 2,
    verbose: bool = False,
) -> List[BenchResult]:
    all_results: List[BenchResult] = []
    for op in ops:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Testing {op} @ {bits}-bit", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)
        res = run_test(op, bits, epochs, binaries=binaries, detail_flag=detail_flag, 
                      repeat=repeat, warmup=warmup, verbose=verbose)
        all_results.extend(res)
    return all_results

def save_results_csv(results: List[BenchResult], path: Union[str, Path]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["op", "bits", "epochs", "metric", "mean_ms", "stdev_ms", 
                       "min_ms", "max_ms", "repeat", "rc", "notes"],
        )
        writer.writeheader()
        for r in results:
            row = {
                "op": r.op,
                "bits": r.bits,
                "epochs": r.epochs,
                "metric": r.metric,
                "mean_ms": f"{r.mean_ms:.3f}",
                "stdev_ms": f"{r.stdev_ms:.3f}",
                "min_ms": f"{r.min_ms:.3f}",
                "max_ms": f"{r.max_ms:.3f}",
                "repeat": r.repeat,
                "rc": r.rc,
                "notes": r.notes,
            }
            writer.writerow(row)

def _print_table(results: List[BenchResult]) -> None:
    if not results:
        print("无结果")
        return
    
    print(f"\n{'='*100}")
    print(f"{'Op':<8} {'Bits':<6} {'Metric':<6} {'Mean(ms)':>10} {'Stdev':>8} {'Min':>8} {'Max':>8} {'N':>4} {'CV%':>6}")
    print(f"{'='*100}")
    
    for r in results:
        mean_str = f"{r.mean_ms:.3f}" if r.mean_ms == r.mean_ms else "nan"
        stdev_str = f"{r.stdev_ms:.3f}" if r.stdev_ms == r.stdev_ms else "nan"
        min_str = f"{r.min_ms:.3f}" if r.min_ms == r.min_ms else "nan"
        max_str = f"{r.max_ms:.3f}" if r.max_ms == r.max_ms else "nan"
        

        if r.mean_ms > 0 and r.stdev_ms == r.stdev_ms:
            cv = (r.stdev_ms / r.mean_ms) * 100
            cv_str = f"{cv:.2f}%"
        else:
            cv_str = "N/A"
        
        print(
            f"{r.op:<8} {r.bits:<6} {r.metric:<6} "
            f"{mean_str:>10} {stdev_str:>8} {min_str:>8} {max_str:>8} "
            f"{r.repeat:>4} {cv_str:>6}"
        )
    
    print(f"{'='*100}\n")

def main():
    parser = argparse.ArgumentParser(
        description="批量运行 Coprocessor-operator 算子并统计性能",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # Benchmark 64-bit SELECT with 20 repetitions and 5 warm-up runs
  python3 bench_ops.py --ops SELECT --bits 64 --repeat 20 --warmup 5 -v
  
  # Use the default operators (ADD, EQ, SELECT, AND)
  python3 bench_ops.py --bits 64 --repeat 10
  
  # Benchmark multiple operators and save the results as CSV
  python3 bench_ops.py --ops GT,EQ,ADD --bits 32,64 --repeat 10 --csv results.csv
  
  # Quick single-run benchmark
  python3 bench_ops.py --ops SELECT --bits 64 --repeat 1 --warmup 0
        """
    )
    
    parser.add_argument("--ops", type=str,
                        default="ADD,GT,EQ,SELECT,OR",
                        help="要测试的算子，逗号分隔（默认: ADD,GT,EQ,SELECT,OR）")
    parser.add_argument("--bits", type=str, required=True,
                        help='bit 宽度: 单值 (如 "128")、逗号列表 (如 "8,16,32") 或 "all"')
    parser.add_argument("--epochs", type=int, default=None, 
                        help="测试轮数（epochs），如未提供则按各程序默认")
    parser.add_argument("--repeat", type=int, default=10,
                        help="每个测试重复次数（默认10次）")
    parser.add_argument("--warmup", type=int, default=2,
                        help="预热次数（默认2次，不计入统计）")
    parser.add_argument("--detail-flag", type=str, default=None, 
                        choices=["--no-epoch-detail", "--show-epochs"],
                        help="是否显示每轮明细")
    parser.add_argument("--csv", type=str, default=None, 
                        help="将结果保存到 CSV 的路径")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="显示每次运行的详细信息")


    for key in DEFAULT_BINARIES.keys():
        parser.add_argument(f"--bin.{key}", dest=f"bin_{key}", type=str, default=None,
                            help=f"覆盖 {key} 的可执行路径")

    args = parser.parse_args()



    override_bins = {}
    for key in DEFAULT_BINARIES.keys():
        v = getattr(args, f"bin_{key}", None)
        if v:
            override_bins[key] = v

    ops = [s.strip().upper() for s in args.ops.split(",") if s.strip()]
    
    results = run_benchmarks(
        ops=ops,
        bits=args.bits,
        epochs=args.epochs,
        binaries=override_bins or None,
        detail_flag=args.detail_flag,
        repeat=args.repeat,
        warmup=args.warmup,
        verbose=args.verbose,
    )

    _print_table(results)

    if args.csv:
        save_results_csv(results, args.csv)
        print(f"✅ CSV 已保存：{args.csv}")

if __name__ == "__main__":
    main()
