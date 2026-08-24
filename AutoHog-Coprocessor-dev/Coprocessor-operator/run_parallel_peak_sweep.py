#!/usr/bin/env python3
"""
parallel peak sweep — 替代 run_parallel_peak_sweep.sh (Python 版)

用法:
    python Coprocessor-operator/run_parallel_peak_sweep.py [选项]

输出在 result/YYYYmmdd_HHMMSS/ 下:
  raw.csv              — 所有运行的原始 CSV
  results.xlsx         — 按 (bit-mode) 分工作表的 Excel (结果1)
  plots/fixed_dag/     — 固定 DAG 的折线图 (结果2)
  plots/heatmaps/      — 热力图 (结果3)

选项:
  --bits BITS           位宽列表，逗号分隔 (默认: 256)
  --modes MODES         输入模式: trivial,std，逗号分隔 (默认: std)
  --ops OPS             算子列表，逗号分隔 (默认: add,sub,and,xor,gt_unsigned,select)
  --variants VARIANTS   变体: serial,parallel，逗号分隔 (默认: parallel)
  --warmup N            预热轮数 (默认: 1)
  --repeat N            计时轮数 (默认: 10)
  --op-threads LIST     op_threads 值，逗号分隔或 auto (默认: auto = 2 的幂到 CORES)
  --dag-parallel LIST   dag_parallel 值，逗号分隔或 auto (默认: auto = 2 的幂到 CORES)
  --max-total-threads N 总线程上限 (默认: CORES * 2)
  --no-serial           跳过 serial baseline
  --executable PATH     可执行文件路径
  --build               cmake 重新构建
  --no-plot             跳过绘图
  --output-dir DIR      输出根目录 (默认: result)
"""
import argparse
import csv
import datetime
import itertools
import multiprocessing
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei',
                                    'WenQuanYi Micro Hei',
                                    'Noto Sans CJK SC',
                                    'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


MARKERS = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', 'x', '+', 'd']
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
          '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5']
SHEET_STYLES: dict[str, dict] = {}

OP_LABELS = {
    'add': 'ADD', 'sub': 'SUB', 'add_unsigned': 'ADD_unsigned',
    'sub_unsigned': 'SUB_unsigned', 'and': 'AND', 'or': 'OR',
    'xor': 'XOR', 'select': 'SELECT',
    'gt': 'GT', 'ge': 'GE', 'lt': 'LT', 'le': 'LE',
    'eq': 'EQ', 'ne': 'NE',
    'gt_unsigned': 'GT_unsigned', 'ge_unsigned': 'GE_unsigned',
    'lt_unsigned': 'LT_unsigned', 'le_unsigned': 'LE_unsigned',
}

XLSX_COLUMNS = [
    'op', 'variant', 'input_mode', 'dag_parallel', 'op_threads',
    'Repeat', 'total_ops', 'total_ms', 'avg_repeat_ms',
    'avg_op_ms', 'throughput_ops_s',
]


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def detect_cores() -> int:
    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        return 1


def gen_levels(max_val: int) -> list[int]:
    levels: list[int] = []
    v = 1
    while v <= max_val:
        levels.append(v)
        v <<= 1
    if levels[-1] != max_val:
        levels.append(max_val)
    return levels


def parse_csv_value(value: str) -> list[int]:
    if value.lower() == 'auto':
        return gen_levels(detect_cores())
    return [int(x.strip()) for x in value.split(',') if x.strip()]


def parse_str_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(',') if x.strip()]


def now_compact() -> str:
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')


def now_timestamp() -> str:
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def make_sheet_styles(bits: list[int], modes: list[str]) -> dict[str, dict]:
    """生成 sheet 名到绘图样式的映射"""
    styles: dict[str, dict] = {}
    for i, (bit, mode) in enumerate(itertools.product(bits, modes)):
        name = f'{bit}-{mode}'
        styles[name] = {
            'label': name,
            'color': COLORS[i % len(COLORS)],
            'marker': MARKERS[i % len(MARKERS)],
            'ls': '-',
        }
    return styles


def sheet_name(bit: int, mode: str) -> str:
    return f'{bit}-{mode}'


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Parallel peak sweep — 运行 bench_parallel_peak 并生成结果')
    p.add_argument('--bits', default='256',
                   help='位宽列表，逗号分隔 (默认: 256)')
    p.add_argument('--modes', default='std',
                   help='输入模式: trivial,std，逗号分隔 (默认: std)')
    p.add_argument('--ops', default='add,sub,and,xor,gt_unsigned,select',
                   help='算子列表，逗号分隔')
    p.add_argument('--variants', default='parallel',
                   help='变体: serial,parallel，逗号分隔 (默认: parallel)')
    p.add_argument('--warmup', type=int, default=1, help='预热轮数 (默认: 1)')
    p.add_argument('--repeat', type=int, default=10, help='计时轮数 (默认: 10)')
    p.add_argument('--op-threads', default='auto',
                   help='op_threads 值，逗号分隔或 auto (默认: auto)')
    p.add_argument('--dag-parallel', default='auto',
                   help='dag_parallel 值，逗号分隔或 auto (默认: auto)')
    p.add_argument('--max-total-threads', type=int, default=0,
                   help='总线程上限 (默认: CORES * 2)')
    p.add_argument('--no-serial', action='store_true', help='跳过 serial baseline')
    p.add_argument('--executable', default='', help='可执行文件路径')
    p.add_argument('--build', action='store_true', help='cmake 重新构建')
    p.add_argument('--no-plot', action='store_true', help='跳过绘图')
    p.add_argument('--output-dir', default='result', help='输出根目录 (默认: result)')
    return p.parse_args(argv)


def resolve_args(args: argparse.Namespace) -> dict:
    cores = detect_cores()
    bits = [int(b) for b in parse_str_list(args.bits)]
    modes = parse_str_list(args.modes)
    ops = parse_str_list(args.ops)
    variants = parse_str_list(args.variants)
    op_threads = parse_csv_value(args.op_threads)
    dag_parallel = parse_csv_value(args.dag_parallel)
    max_total = args.max_total_threads if args.max_total_threads > 0 else cores * 2

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.normpath(os.path.join(script_dir, '..'))
    if args.executable:
        exe = os.path.abspath(args.executable)
    else:
        exe = os.path.join(repo_root, 'build', 'Coprocessor-operator', 'bench_parallel_peak')
        exe = os.path.abspath(exe)

    return {
        'cores': cores,
        'bits': bits,
        'modes': modes,
        'ops': ops,
        'variants': variants,
        'op_threads': op_threads,
        'dag_parallel': dag_parallel,
        'max_total_threads': max_total,
        'warmup': args.warmup,
        'repeat': args.repeat,
        'run_serial': not args.no_serial,
        'exe': exe,
        'build': args.build,
        'no_plot': args.no_plot,
        'output_root': args.output_dir,
        'repo_root': repo_root,
    }


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def do_build(cfg: dict) -> None:
    repo = cfg['repo_root']
    build_dir = os.path.join(repo, 'build')
    cores = cfg['cores']
    print(f'[build] cmake -S {repo} -B {build_dir} -DCMAKE_BUILD_TYPE=Release')
    subprocess.run(['cmake', '-S', repo, '-B', build_dir,
                    '-DCMAKE_BUILD_TYPE=Release'], check=True)
    print(f'[build] cmake --build {build_dir} --target bench_parallel_peak -j{cores}')
    subprocess.run(['cmake', '--build', build_dir,
                    '--target', 'bench_parallel_peak', '-j', str(cores)], check=True)
    print('[build] 完成')


def run_one(cfg: dict, run_log_dir: str,
            width: int, op: str, variant: str, input_mode: str,
            op_threads: int, dag_parallel: int) -> dict | None:
    """运行一次 bench_parallel_peak，返回解析后的字典"""
    exe = cfg['exe']
    log_file = os.path.join(
        run_log_dir,
        f'run_{input_mode}_w{width}_{op}_{variant}_ot{op_threads}_dp{dag_parallel}.log')

    cmd = [
        exe, '--quiet',
        '--width', str(width),
        '--op', op,
        '--variant', variant,
        '--op-threads', str(op_threads),
        '--dag-parallel', str(dag_parallel),
        '--warmup', str(cfg['warmup']),
        '--repeat', str(cfg['repeat']),
        '--input-mode', input_mode,
    ]

    env = os.environ.copy()
    env['OMP_DYNAMIC'] = 'false'
    env['OMP_PROC_BIND'] = 'spread'
    env['OMP_PLACES'] = 'cores'
    env['OMP_NUM_THREADS'] = str(op_threads)

    print(f'  [{input_mode} w{width} {op} {variant} '
          f'ot{op_threads} dp{dag_parallel}]', end=' ', flush=True)

    try:
        result = subprocess.run(cmd, env=env, capture_output=True,
                                text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        print('TIMEOUT')
        return None

    if result.returncode != 0:
        print(f'FAILED (rc={result.returncode})')
        with open(log_file, 'w') as f:
            f.write('STDOUT:\n' + result.stdout + '\nSTDERR:\n' + result.stderr)
        return None

    line = result.stdout.strip()
    if not line:
        print('EMPTY OUTPUT')
        return None


    # timestamp,width,op,variant,input_mode,op_threads,dag_parallel,
    # warmup,repeat,total_ops,total_ms,avg_repeat_ms,avg_op_ms,
    # throughput_ops_s,checksum,hw_threads,omp_max_threads
    parts = line.split(',')
    if len(parts) < 14:
        print(f'PARSE ERROR (got {len(parts)} fields)')
        with open(log_file, 'w') as f:
            f.write(result.stdout)
        return None

    try:
        rec = {
            'width': int(parts[1]),
            'op': parts[2],
            'variant': parts[3],
            'input_mode': parts[4],
            'op_threads': int(parts[5]),
            'dag_parallel': int(parts[6]),
            'Repeat': int(parts[8]),
            'total_ops': int(parts[9]),
            'total_ms': float(parts[10]),
            'avg_repeat_ms': float(parts[11]),
            'avg_op_ms': float(parts[12]),
            'throughput_ops_s': float(parts[13]),
        }
    except (ValueError, IndexError) as e:
        print(f'PARSE ERROR: {e}')
        with open(log_file, 'w') as f:
            f.write(result.stdout + '\n' + result.stderr)
        return None

    print(f'{rec["throughput_ops_s"]:.2f} ops/s')
    with open(log_file, 'w') as f:
        f.write(result.stdout)
        if result.stderr:
            f.write('\nSTDERR:\n' + result.stderr)
    return rec


def run_sweep(cfg: dict, out_dir: str) -> list[dict]:
    """执行完整 sweep"""
    run_log_dir = os.path.join(out_dir, 'run_logs')
    os.makedirs(run_log_dir, exist_ok=True)

    results: list[dict] = []


    total = 0
    for width, mode, op in itertools.product(cfg['bits'], cfg['modes'], cfg['ops']):
        if cfg['run_serial']:
            total += 1
        for variant in cfg['variants']:
            if variant == 'serial':
                continue
            for dp, ot in itertools.product(cfg['dag_parallel'], cfg['op_threads']):
                if ot * dp > cfg['max_total_threads']:
                    continue
                total += 1

    print(f'[info] cores={cfg["cores"]}')
    print(f'[info] bits={cfg["bits"]}')
    print(f'[info] modes={cfg["modes"]}')
    print(f'[info] ops={cfg["ops"]}')
    print(f'[info] variants={cfg["variants"]}')
    print(f'[info] op_threads={cfg["op_threads"]}')
    print(f'[info] dag_parallel={cfg["dag_parallel"]}')
    print(f'[info] max_total_threads={cfg["max_total_threads"]}')
    print(f'[info] warmup={cfg["warmup"]} repeat={cfg["repeat"]}')
    print(f'[info] total_runs={total}')
    print()

    count = 0
    for width, mode, op in itertools.product(cfg['bits'], cfg['modes'], cfg['ops']):
        # serial baseline
        if cfg['run_serial']:
            count += 1
            print(f'[{count}/{total}] ', end='')
            rec = run_one(cfg, run_log_dir, width, op, 'serial', mode, 1, 1)
            if rec:
                results.append(rec)

        for variant in cfg['variants']:
            if variant == 'serial':
                continue
            for dp, ot in itertools.product(cfg['dag_parallel'], cfg['op_threads']):
                if ot * dp > cfg['max_total_threads']:
                    continue
                count += 1
                print(f'[{count}/{total}] ', end='')
                rec = run_one(cfg, run_log_dir, width, op, variant, mode, ot, dp)
                if rec:
                    results.append(rec)

    print(f'\n[done] {len(results)}/{total} runs succeeded')
    return results


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def write_xlsx(results: list[dict], out_dir: str,
               bits: list[int], modes: list[str]) -> str:
    """
    生成 results.xlsx，按 (bit-mode) 分工作表。
    表头: op, variant, input_mode, dag_parallel, op_threads,
          Repeat, total_ops, total_ms, avg_repeat_ms, avg_op_ms, throughput_ops_s
    """
    from openpyxl import Workbook
    wb = Workbook()
    wb.remove(wb.active)


    groups: dict[str, list[dict]] = defaultdict(list)
    for rec in results:
        key = sheet_name(rec['width'], rec['input_mode'])
        groups[key].append(rec)

    for bit in bits:
        for mode in modes:
            sname = sheet_name(bit, mode)
            ws = wb.create_sheet(title=sname)
            ws.append(XLSX_COLUMNS)

            sorted_recs = sorted(groups.get(sname, []),
                                 key=lambda r: (r['dag_parallel'], r['op_threads']))
            for rec in sorted_recs:
                row = [rec.get(col, '') for col in XLSX_COLUMNS]
                ws.append(row)

    xlsx_path = os.path.join(out_dir, 'results.xlsx')
    wb.save(xlsx_path)
    print(f'[xlsx] 已保存: {xlsx_path}')
    print(f'      工作表: {[s.title for s in wb.worksheets]}')
    return xlsx_path


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def plot_fixed_dag(results: list[dict], out_dir: str,
                   bits: list[int], modes: list[str],
                   dag_value: int = 1) -> None:
    """
    结果2: 固定 DAG (dag_parallel = dag_value),
    X 轴 = OMP 线程 (op_threads),
    不同 (bit-mode) 组合作为不同线条.
    每个算子一张图 + 一张汇总大图.
    """
    styles = make_sheet_styles(bits, modes)
    plots_dir = os.path.join(out_dir, 'plots', 'fixed_dag')
    os.makedirs(plots_dir, exist_ok=True)

    all_ops = sorted(set(r['op'] for r in results))


    for op_name in all_ops:
        fig, ax = plt.subplots(figsize=(9, 6))
        op_label = OP_LABELS.get(op_name, op_name)
        fig.suptitle(f'{op_label}  —  Throughput vs OMP  (DAG={dag_value})',
                     fontsize=14, fontweight='bold')

        has_data = False
        for bit in bits:
            for mode in modes:
                sname = sheet_name(bit, mode)
                style = styles.get(sname)
                if not style:
                    continue


                recs = [r for r in results
                        if r['width'] == bit
                        and r['input_mode'] == mode
                        and r['op'] == op_name
                        and r['dag_parallel'] == dag_value
                        and r['variant'] != 'serial']
                if not recs:
                    continue

                recs.sort(key=lambda x: x['op_threads'])
                omp_vals = [r['op_threads'] for r in recs]
                tp_vals = [r['throughput_ops_s'] for r in recs]

                ax.plot(omp_vals, tp_vals,
                        color=style['color'], marker=style['marker'],
                        linestyle=style['ls'], linewidth=2.2,
                        label=style['label'], markersize=7,
                        markeredgewidth=1, markeredgecolor='black')
                has_data = True

        if not has_data:
            plt.close(fig)
            continue

        ax.set_xlabel('OMP Threads (OpenMP Parallelism)', fontsize=12)
        ax.set_ylabel('Throughput (ops/s)', fontsize=12)
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xscale('log', base=2)

        plt.tight_layout()
        path = os.path.join(plots_dir, f'{op_name}_fixed_dag{dag_value}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  [固定DAG] {op_label} -> {path}')


    n_ops = len(all_ops)
    if n_ops == 0:
        return
    n_cols = min(4, n_ops)
    n_rows = (n_ops + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 4.5 * n_rows),
                             squeeze=False)
    fig.suptitle(f'Fixed DAG={dag_value}: Throughput vs OMP Threads',
                 fontsize=14, fontweight='bold')

    for idx, op_name in enumerate(all_ops):
        r, c = divmod(idx, n_cols)
        ax = axes[r][c]

        for bit in bits:
            for mode in modes:
                sname = sheet_name(bit, mode)
                style = styles.get(sname)
                if not style:
                    continue
                recs = [r for r in results
                        if r['width'] == bit
                        and r['input_mode'] == mode
                        and r['op'] == op_name
                        and r['dag_parallel'] == dag_value
                        and r['variant'] != 'serial']
                if not recs:
                    continue
                recs.sort(key=lambda x: x['op_threads'])
                omp_vals = [r['op_threads'] for r in recs]
                tp_vals = [r['throughput_ops_s'] for r in recs]
                ax.plot(omp_vals, tp_vals,
                        color=style['color'], marker=style['marker'],
                        linestyle=style['ls'], linewidth=1.5,
                        label=style['label'], markersize=4)

        ax.set_xlabel('OMP')
        ax.set_ylabel('Throughput')
        ax.set_title(f'{OP_LABELS.get(op_name, op_name)}',
                     fontsize=10, fontweight='bold')
        ax.legend(fontsize=6, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xscale('log', base=2)

    for idx in range(n_ops, n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r][c].set_visible(False)

    plt.tight_layout()
    overview_path = os.path.join(plots_dir, f'overview_fixed_dag{dag_value}.png')
    fig.savefig(overview_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  [固定DAG汇总] -> {overview_path}')


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def plot_heatmaps(results: list[dict], out_dir: str,
                  bits: list[int], modes: list[str]) -> None:
    """
    结果3: 热力图.
    每个 (bit, mode, op) 组合一张热力图:
      X 轴 = dag_parallel, Y 轴 = op_threads, 颜色 = throughput
    """
    heat_dir = os.path.join(out_dir, 'plots', 'heatmaps')
    os.makedirs(heat_dir, exist_ok=True)

    all_ops = sorted(set(r['op'] for r in results))

    for bit in bits:
        for mode in modes:

            ops_present = [
                op for op in all_ops
                if any(r['width'] == bit and r['input_mode'] == mode
                       and r['op'] == op for r in results)
            ]
            if not ops_present:
                continue

            n_ops = len(ops_present)
            n_cols = min(4, n_ops)
            n_rows = (n_ops + n_cols - 1) // n_cols

            fig, axes = plt.subplots(n_rows, n_cols,
                                     figsize=(5 * n_cols, 4 * n_rows),
                                     squeeze=False)
            fig.suptitle(f'Throughput Heatmap  ({bit}-{mode})',
                         fontsize=14, fontweight='bold')

            for idx, op_name in enumerate(ops_present):
                r, c = divmod(idx, n_cols)
                ax = axes[r][c]


                recs = [rec for rec in results
                        if rec['width'] == bit
                        and rec['input_mode'] == mode
                        and rec['op'] == op_name
                        and rec['variant'] != 'serial']
                if not recs:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                            transform=ax.transAxes)
                    ax.set_title(f'{OP_LABELS.get(op_name, op_name)}')
                    continue

                threads_set = sorted(set(r['op_threads'] for r in recs))
                dags_set = sorted(set(r['dag_parallel'] for r in recs))

                if len(threads_set) < 2 or len(dags_set) < 2:
                    ax.text(0.5, 0.5,
                            f'Insufficient data\n({len(threads_set)} OMP × '
                            f'{len(dags_set)} DAG)',
                            ha='center', va='center', transform=ax.transAxes,
                            fontsize=10, style='italic', color='gray')
                    ax.set_title(f'{OP_LABELS.get(op_name, op_name)}')
                    ax.set_xlabel('DAG Parallelism')
                    ax.set_ylabel('OMP Threads')
                    continue

                tp_map: dict[tuple[int, int], float] = {}
                for rec in recs:
                    tp_map[(rec['op_threads'], rec['dag_parallel'])] = \
                        rec['throughput_ops_s']

                matrix = np.zeros((len(threads_set), len(dags_set)))
                for i, t in enumerate(threads_set):
                    for j, d in enumerate(dags_set):
                        matrix[i, j] = tp_map.get((t, d), 0.0)

                im = ax.imshow(matrix, aspect='auto', cmap='viridis')
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                             label='Throughput (ops/s)')

                ax.set_xticks(range(len(dags_set)))
                ax.set_xticklabels([str(d) for d in dags_set], fontsize=8)
                ax.set_yticks(range(len(threads_set)))
                ax.set_yticklabels([str(t) for t in threads_set], fontsize=8)
                ax.set_xlabel('DAG Parallelism')
                ax.set_ylabel('OMP Threads')
                ax.set_title(f'{OP_LABELS.get(op_name, op_name)}')


                for i in range(len(threads_set)):
                    for j in range(len(dags_set)):
                        val = matrix[i, j]
                        if val > 0:
                            text_color = 'white' if val < matrix.max() * 0.7 else 'black'
                            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                                    fontsize=6, color=text_color)

            for idx in range(n_ops, n_rows * n_cols):
                r, c = divmod(idx, n_cols)
                axes[r][c].set_visible(False)

            plt.tight_layout()
            out_path = os.path.join(heat_dir, f'heatmap_{bit}-{mode}.png')
            fig.savefig(out_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f'  [热力图] {bit}-{mode} -> {out_path}')


    for bit in bits:
        for mode in modes:
            for op_name in all_ops:
                recs = [rec for rec in results
                        if rec['width'] == bit
                        and rec['input_mode'] == mode
                        and rec['op'] == op_name
                        and rec['variant'] != 'serial']
                if not recs:
                    continue

                threads_set = sorted(set(r['op_threads'] for r in recs))
                dags_set = sorted(set(r['dag_parallel'] for r in recs))
                if len(threads_set) < 2 or len(dags_set) < 2:
                    continue

                tp_map: dict[tuple[int, int], float] = {}
                for rec in recs:
                    tp_map[(rec['op_threads'], rec['dag_parallel'])] = \
                        rec['throughput_ops_s']

                matrix = np.zeros((len(threads_set), len(dags_set)))
                for i, t in enumerate(threads_set):
                    for j, d in enumerate(dags_set):
                        matrix[i, j] = tp_map.get((t, d), 0.0)

                fig, ax = plt.subplots(figsize=(8, 6))
                op_label = OP_LABELS.get(op_name, op_name)
                fig.suptitle(
                    f'{op_label}  —  Throughput Heatmap  ({bit}-{mode})',
                    fontsize=14, fontweight='bold')

                im = ax.imshow(matrix, aspect='auto', cmap='viridis')
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                             label='Throughput (ops/s)')

                ax.set_xticks(range(len(dags_set)))
                ax.set_xticklabels([str(d) for d in dags_set], fontsize=10)
                ax.set_yticks(range(len(threads_set)))
                ax.set_yticklabels([str(t) for t in threads_set], fontsize=10)
                ax.set_xlabel('DAG Parallelism', fontsize=12)
                ax.set_ylabel('OMP Threads', fontsize=12)

                for i in range(len(threads_set)):
                    for j in range(len(dags_set)):
                        val = matrix[i, j]
                        if val > 0:
                            text_color = 'white' if val < matrix.max() * 0.7 else 'black'
                            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                                    fontsize=8, color=text_color)

                plt.tight_layout()
                ind_dir = os.path.join(heat_dir, 'individual')
                os.makedirs(ind_dir, exist_ok=True)
                path = os.path.join(ind_dir,
                                    f'heatmap_{bit}-{mode}_{op_name}.png')
                fig.savefig(path, dpi=150, bbox_inches='tight')
                plt.close(fig)


# ═══════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════

def main() -> None:
    args = parse_args()
    cfg = resolve_args(args)

    print('═' * 60)
    print('  AutoHog parallel peak sweep')
    print(f'  {now_timestamp()}')
    print(f'  cores={cfg["cores"]}')
    print('═' * 60)
    print()

    if cfg['build']:
        do_build(cfg)

    if not os.path.isfile(cfg['exe']):
        print(f'[error] 可执行文件不存在: {cfg["exe"]}')
        print('        请先构建: cmake --build build --target bench_parallel_peak')
        sys.exit(1)

    ts = now_compact()
    out_dir = os.path.join(cfg['output_root'], ts)
    os.makedirs(out_dir, exist_ok=True)
    print(f'[output] {out_dir}')
    print()


    results = run_sweep(cfg, out_dir)

    if not results:
        print('[error] 所有 run 都失败了')
        sys.exit(1)


    raw_csv = os.path.join(out_dir, 'raw.csv')
    fieldnames = [
        'width', 'op', 'variant', 'input_mode', 'op_threads', 'dag_parallel',
        'Repeat', 'total_ops', 'total_ms', 'avg_repeat_ms',
        'avg_op_ms', 'throughput_ops_s',
    ]
    with open(raw_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    print(f'[csv] {raw_csv}  ({len(results)} rows)')


    print('[xlsx] 生成中...')
    write_xlsx(results, out_dir, cfg['bits'], cfg['modes'])
    print()


    if not cfg['no_plot']:
        print('[plot] 固定 DAG 折线图 (结果2)...')

        dag_values = sorted(set(r['dag_parallel'] for r in results))
        for dag_val in dag_values:
            print(f'  DAG={dag_val} ...')
            plot_fixed_dag(results, out_dir, cfg['bits'], cfg['modes'],
                           dag_value=dag_val)

        print('[plot] 热力图 (结果3)...')
        plot_heatmaps(results, out_dir, cfg['bits'], cfg['modes'])

        print('[plot] 完成')
        print()

    print(f'[done] 所有结果在 {out_dir}/')


if __name__ == '__main__':
    main()
