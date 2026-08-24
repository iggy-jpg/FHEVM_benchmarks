#!/bin/bash

#



#




#   AUTOHOG_PERF_OPS="add_u sub_u select mul_plain_low div_plain_q" ./run_perf_all.sh
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/build/Coprocessor-operator/test_PERF_ALL"
WIDTHS=${1:-"64 128"}
EPOCHS=${2:-2}
COOL=${3:-8}
OUT="$ROOT/results/perf_all_$(date +%Y%m%d_%H%M%S).txt"

if [ ! -x "$BIN" ]; then
  echo "未找到 $BIN；请先构建： cmake --build build --target test_PERF_ALL -j" >&2
  exit 1
fi
mkdir -p "$ROOT/results"


OPS=${AUTOHOG_PERF_OPS:-$($BIN 0 0 list)}
echo "算子: $OPS" | tee "$OUT"
echo "线程: ${OMP_NUM_THREADS:-auto}  轮数: $EPOCHS  冷却: ${COOL}s" | tee -a "$OUT"

for W in $WIDTHS; do
  for OP in $OPS; do

    "$BIN" "$W" "$EPOCHS" "$OP" 2>&1 | grep -E "^PERF|FAIL" | tee -a "$OUT"
    sleep "$COOL"
  done
done

echo "" | tee -a "$OUT"
echo "===== 汇总（毫秒，只含并行变体） =====" | tee -a "$OUT"
python3 - "$OUT" <<'PYEOF' | tee -a "$OUT"
import re, sys, collections
rows = collections.defaultdict(dict)
cls = {}
for line in open(sys.argv[1]):
    m = re.match(r'PERF op=(\S+) W=(\d+) ms=([\d.]+) pbs_class=(\S+) correct=(\d+)/(\d+)', line)
    if not m: continue
    op, W, ms, c, ok, tot = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4), int(m.group(5)), int(m.group(6))
    rows[op][W] = (ms, ok == tot)
    cls[op] = c
widths = sorted({w for v in rows.values() for w in v})
order = {'ct': 0, 'plain': 1, 'free': 2}
hdr = f"{'算子':<16}{'类别':<8}" + "".join(f"{'W='+str(w):>14}" for w in widths) + "  正确"
print(hdr); print('-' * len(hdr))
for op in sorted(rows, key=lambda o: (order.get(cls[o], 9), o)):
    line = f"{op:<16}{cls[op]:<8}"
    allok = True
    for w in widths:
        if w in rows[op]:
            ms, ok = rows[op][w]
            line += f"{ms:>13.1f}ms"
            allok = allok and ok
        else:
            line += f"{'-':>14}"
    print(line + ("  ✓" if allok else "  ✗"))
PYEOF
echo "" | tee -a "$OUT"
echo "结果已保存: $OUT"
