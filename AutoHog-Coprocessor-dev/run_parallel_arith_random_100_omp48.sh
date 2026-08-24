#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN="$ROOT/build/Coprocessor-operator/test_PARALLEL_ARITH_RANDOM"
REPEAT=${REPEAT:-100}
THREADS=${THREADS:-48}
WIDTHS=${WIDTHS:-4,8,16,32,64,128,256}
OPS=${OPS:-add,gt,mul,div}
SEED=${SEED:-0x20260726}

if [[ ! -x "$BIN" ]]; then
  echo "missing executable: $BIN" >&2
  exit 1
fi

TS=$(date +%Y%m%d_%H%M%S)
OUT_DIR="$ROOT/results/parallel_arith_random_r${REPEAT}_omp${THREADS}_${TS}"
mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/run.log"
RAW="$OUT_DIR/raw.csv"
SUMMARY="$OUT_DIR/summary.csv"
STATUS="$OUT_DIR/run.status"
PIDFILE="$OUT_DIR/run.pid"

printf 'state=starting\nstarted_at=%s\n' "$(date '+%F %T')" > "$STATUS"

nohup bash -lc '
  set +e
  status=$1 threads=$2 bin=$3 widths=$4 ops=$5 repeat=$6 seed=$7
  raw=$8 summary=$9 out_dir=${10}
  started=$(date "+%F %T")
  printf "state=running\nstarted_at=%s\n" "$started" > "$status"
  printf "width,op,sample,threads,time_ms,correct,rss_kib,peak_rss_kib,process_threads,finished_at,status\n" > "$raw"
  printf "width,op,variant,input_mode,threads,repeat,correct,total,total_ms,avg_ms,min_ms,max_ms,peak_rss_kib\n" > "$summary"
  rc=0
  IFS=, read -r -a width_list <<< "$widths"
  for width in "${width_list[@]}"; do
    printf "state=running\nstarted_at=%s\ncurrent_width=%s\n" "$started" "$width" > "$status"
    width_raw="$out_dir/raw_w${width}.csv"
    width_summary="$out_dir/summary_w${width}.csv"
    resource="$out_dir/resource_w${width}.txt"
    echo "[WIDTH_START] $(date "+%F %T") width=$width"
    OMP_NUM_THREADS="$threads" OMP_THREAD_LIMIT="$threads" OMP_DYNAMIC=FALSE \
      /usr/bin/time -v -o "$resource" stdbuf -oL -eL "$bin" \
        --widths "$width" --ops "$ops" --repeat "$repeat" --threads "$threads" \
        --progress-every 1 --seed "$seed" --raw-csv "$width_raw" \
        --summary-csv "$width_summary"
    width_rc=$?
    if [[ -s "$width_raw" ]]; then tail -n +2 "$width_raw" >> "$raw"; fi
    if [[ -s "$width_summary" ]]; then tail -n +2 "$width_summary" >> "$summary"; fi
    echo "[WIDTH_DONE] $(date "+%F %T") width=$width exit_code=$width_rc"
    if [[ $width_rc -ne 0 ]]; then rc=$width_rc; break; fi
  done
  if [[ $rc -eq 0 ]]; then state=completed; else state=failed; fi
  printf "state=%s\nfinished_at=%s\nexit_code=%s\n" \
    "$state" "$(date "+%F %T")" "$rc" > "$1"
  exit "$rc"
' bash "$STATUS" "$THREADS" "$BIN" "$WIDTHS" "$OPS" "$REPEAT" \
  "$SEED" "$RAW" "$SUMMARY" "$OUT_DIR" > "$LOG" 2>&1 &

PID=$!
printf '%s\n' "$PID" > "$PIDFILE"
printf 'pid=%s\noutput_dir=%s\nlog=%s\nstatus=%s\nraw_csv=%s\nsummary_csv=%s\n' \
  "$PID" "$OUT_DIR" "$LOG" "$STATUS" "$RAW" "$SUMMARY"
