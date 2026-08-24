#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mem_watch_add.py (单位：MB)
依次运行不同位宽（8,16,32,64,128,256）的 ./Coprocessor-operator/test_ADD <W> <reps>
采样进程RSS内存(单位MB)并写入CSV，最终生成曲线图 PNG。
输出路径自动为 mem_logs/<执行文件名>/
"""

import argparse
import os
import time
import subprocess
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt


def read_rss_mb_linux(pid: int) -> float:
    """从 /proc/<pid>/status 读取 VmRSS (kB)，转换为 MB（MiB）。"""
    with open(f"/proc/{pid}/status", "r") as f:
        for line in f:
            if line.startswith("VmRSS:"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    rss_kb = int(parts[1])
                    return rss_kb / 1024.0
    return 0.0


def run_and_log_one(exec_path: str, bitwidth: int, reps: int, interval: float, outdir: Path):
    """运行单个位宽测试并记录RSS(MB)。"""
    csv_path = outdir / f"mem_W{bitwidth}.csv"
    log_path = outdir / f"child_W{bitwidth}.log"

    cmd = [exec_path, str(bitwidth), str(reps)]
    print(f"[INFO] Running: {' '.join(cmd)}")
    start_wall = time.time()

    with open(log_path, "w") as child_log:
        proc = subprocess.Popen(cmd, stdout=child_log, stderr=subprocess.STDOUT)
        times, rss_vals_mb = [], []

        try:
            while True:
                ret = proc.poll()
                elapsed = time.time() - start_wall
                try:
                    rss_mb = read_rss_mb_linux(proc.pid)
                except (FileNotFoundError, PermissionError):
                    rss_mb = None
                if rss_mb is not None:
                    times.append(elapsed)
                    rss_vals_mb.append(rss_mb)
                if ret is not None:
                    break
                time.sleep(interval)
        except KeyboardInterrupt:
            print("[WARN] Ctrl+C detected, terminating child...")
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()

    if not times:
        times, rss_vals_mb = [0.0], [0.0]

    with open(csv_path, "w") as f:
        f.write("time_s,rss_mb\n")
        for t, r in zip(times, rss_vals_mb):
            f.write(f"{t:.3f},{r:.3f}\n")

    print(f"[OK] W={bitwidth} done. CSV => {csv_path.name}, LOG => {log_path.name}")
    return times, rss_vals_mb


def main():
    parser = argparse.ArgumentParser(description="Run test_ADD with multiple bitwidths and log memory usage (MB).")
    parser.add_argument("--exec", default="./build/Coprocessor-operator/test_ADD", help="可执行文件路径")
    parser.add_argument("--bitwidths", default="8,16,32,64,128,256", help="位宽列表，逗号分隔")
    parser.add_argument("--reps", type=int, default=10, help="每个位宽的num_test/循环次数")
    parser.add_argument("--interval", type=float, default=0.5, help="采样间隔（秒）")
    parser.add_argument("--outdir", default="mem_logs", help="基础输出目录")
    args = parser.parse_args()

    exec_path = Path(args.exec)
    exec_name = exec_path.name


    base_outdir = Path(args.outdir)
    outdir = base_outdir / exec_name
    outdir.mkdir(parents=True, exist_ok=True)

    bitwidths = [int(x) for x in args.bitwidths.split(",") if x.strip()]
    reps = args.reps
    interval = args.interval


    meta_path = outdir / "RUN_INFO.txt"
    with open(meta_path, "a") as mf:
        mf.write(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] exec={exec_path}, "
                 f"bitwidths={bitwidths}, reps={reps}, interval={interval}, unit=MB\n")

    all_curves = {}
    try:
        for w in bitwidths:
            t, r = run_and_log_one(str(exec_path), w, reps, interval, outdir)
            all_curves[w] = (t, r)
    except KeyboardInterrupt:
        print("[WARN] Interrupted by user. Plotting what we have...")

    if all_curves:
        plt.figure(figsize=(8, 5))
        for w in sorted(all_curves.keys()):
            t, r = all_curves[w]
            plt.plot(t, r, label=f"W={w}")
        plt.xlabel("Time (s)")
        plt.ylabel("RSS (MB)")
        plt.title(f"Memory (RSS) over Time for {exec_name}")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        png_path = outdir / "mem_curves.png"
        plt.tight_layout()
        plt.savefig(png_path, dpi=150)
        print(f"[OK] Plot saved => {png_path}")

    print(f"[DONE] All tests completed. Results stored in: {outdir.resolve()}")


if __name__ == "__main__":
    main()