#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXECUTABLE=${EXECUTABLE:-"$ROOT_DIR/build-consensus/Coprocessor-operator/test_consensus"}
CONCURRENCY_LEVELS=(1 2 4 8 16)
REGENERATE_KEYS=false
COMPUTATIONS_FILE=""
DISABLE_PLAIN_FOLD=false
DISABLE_LOW4_SANITIZE=false
SHARED_KEY=false

usage() {
    echo "Usage: $0 [--regenerate-keys] [--executable PATH] [--levels '1 2 4 8 16'] [--computations FILE] [--disable-plain-fold] [--disable-low4-sanitize] [--shared-key]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --regenerate-keys)
            REGENERATE_KEYS=true
            shift
            ;;
        --executable)
            [[ $# -ge 2 ]] || { usage >&2; exit 2; }
            EXECUTABLE=$2
            shift 2
            ;;
        --levels)
            [[ $# -ge 2 ]] || { usage >&2; exit 2; }
            read -r -a CONCURRENCY_LEVELS <<< "$2"
            [[ ${#CONCURRENCY_LEVELS[@]} -gt 0 ]] || { usage >&2; exit 2; }
            shift 2
            ;;
        --computations)
            [[ $# -ge 2 ]] || { usage >&2; exit 2; }
            COMPUTATIONS_FILE=$2
            shift 2
            ;;
        --disable-plain-fold)
            DISABLE_PLAIN_FOLD=true
            shift
            ;;
        --disable-low4-sanitize)
            DISABLE_LOW4_SANITIZE=true
            shift
            ;;
        --shared-key)
            SHARED_KEY=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -x "$EXECUTABLE" ]] || {
    echo "Executable not found: $EXECUTABLE" >&2
    exit 1
}

cd "$ROOT_DIR"
CONSENSUS_WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/autohog-consensus.XXXXXX")
cleanup_work_dir() {
    rm -rf -- "$CONSENSUS_WORK_DIR"
}
trap cleanup_work_dir EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
export AUTOHOG_CONSENSUS_WORKDIR="$CONSENSUS_WORK_DIR"
export AUTOHOG_FFTW_WISDOM_FILE="$CONSENSUS_WORK_DIR/fftw.wisdom"
export AUTOHOG_CONSENSUS_TOTAL_THREADS=${AUTOHOG_CONSENSUS_TOTAL_THREADS:-16}
export OMP_NUM_THREADS=$AUTOHOG_CONSENSUS_TOTAL_THREADS
export OMP_THREAD_LIMIT=$AUTOHOG_CONSENSUS_TOTAL_THREADS
export OMP_DYNAMIC=FALSE
if [[ -n "$COMPUTATIONS_FILE" ]]; then
    [[ -f "$COMPUTATIONS_FILE" ]] || {
        echo "Computations file not found: $COMPUTATIONS_FILE" >&2
        exit 1
    }
    export AUTOHOG_CONSENSUS_COMPUTATIONS="$COMPUTATIONS_FILE"
fi
if [[ "$DISABLE_PLAIN_FOLD" == true ]]; then
    export AUTOHOG_DISABLE_PLAIN_FOLD=1
fi
if [[ "$DISABLE_LOW4_SANITIZE" == true ]]; then
    export AUTOHOG_DISABLE_LOW4_SANITIZE=1
fi
if [[ "$SHARED_KEY" == true ]]; then
    export AUTOHOG_CONSENSUS_SHARED_KEY=1
fi

RUN_ID=$(date +%Y%m%d_%H%M%S)
RESULT_DIR="$ROOT_DIR/results/consensus_raw_${RUN_ID}"
BASELINE_DIR="$RESULT_DIR/baseline_concurrency_1"
INPUT_DIR="$RESULT_DIR/fixed_inputs"
COMBINED_CSV="$RESULT_DIR/summary.csv"
OPERATION_CSV="$RESULT_DIR/operation_summary.csv"
mkdir -p "$RESULT_DIR" "$BASELINE_DIR" "$INPUT_DIR"

if [[ "$REGENERATE_KEYS" == true ]]; then
    rm -f client_secret.bin server_eval.bin server_annihilate.bin server_pk.bin server_pk_v2.bin
fi

export AUTOHOG_CONSENSUS_BASELINE="$BASELINE_DIR"
export AUTOHOG_CONSENSUS_INPUTS="$INPUT_DIR"

echo "run_id=$RUN_ID" > "$RESULT_DIR/run.status"
echo "status=running" >> "$RESULT_DIR/run.status"
echo "started_at=$(date '+%F %T %z')" >> "$RESULT_DIR/run.status"
echo "executable=$EXECUTABLE" >> "$RESULT_DIR/run.status"
echo "omp_num_threads=$OMP_NUM_THREADS" >> "$RESULT_DIR/run.status"
echo "consensus_total_threads=$AUTOHOG_CONSENSUS_TOTAL_THREADS" >> "$RESULT_DIR/run.status"
echo "computations=${AUTOHOG_CONSENSUS_COMPUTATIONS:-data/computations.txt}" >> "$RESULT_DIR/run.status"
echo "plain_fold_enabled=$([[ "$DISABLE_PLAIN_FOLD" == true ]] && echo 0 || echo 1)" >> "$RESULT_DIR/run.status"
echo "low4_sanitize_enabled=$([[ "$DISABLE_LOW4_SANITIZE" == true ]] && echo 0 || echo 1)" >> "$RESULT_DIR/run.status"
echo "shared_key=$([[ "$SHARED_KEY" == true ]] && echo 1 || echo 0)" >> "$RESULT_DIR/run.status"
echo "temporary_work_directory=$AUTOHOG_CONSENSUS_WORKDIR" >> "$RESULT_DIR/run.status"
echo "temporary_fftw_wisdom=$AUTOHOG_FFTW_WISDOM_FILE" >> "$RESULT_DIR/run.status"

for concurrency in "${CONCURRENCY_LEVELS[@]}"; do
    LEVEL_LOG="$RESULT_DIR/concurrency_${concurrency}.log"
    LEVEL_CSV="$RESULT_DIR/concurrency_${concurrency}.csv"
    export AUTOHOG_CONSENSUS_SUMMARY="$LEVEL_CSV"

    echo "[$(date '+%F %T')] starting external concurrency=$concurrency"
    set +e
    if command -v stdbuf >/dev/null 2>&1; then
        stdbuf -oL -eL "$EXECUTABLE" "$concurrency" 2>&1 | tee "$LEVEL_LOG"
        exit_code=${PIPESTATUS[0]}
    else
        "$EXECUTABLE" "$concurrency" 2>&1 | tee "$LEVEL_LOG"
        exit_code=${PIPESTATUS[0]}
    fi
    set -e

    if [[ $exit_code -ne 0 ]]; then
        {
            echo "status=failed"
            echo "failed_concurrency=$concurrency"
            echo "exit_code=$exit_code"
            echo "finished_at=$(date '+%F %T %z')"
        } >> "$RESULT_DIR/run.status"
        echo "Consensus test failed at concurrency=$concurrency (exit=$exit_code)" >&2
        exit "$exit_code"
    fi
    echo "[$(date '+%F %T')] completed external concurrency=$concurrency"
done

if [[ -f "$AUTOHOG_FFTW_WISDOM_FILE" ]]; then
    cp "$AUTOHOG_FFTW_WISDOM_FILE" "$RESULT_DIR/fftw.wisdom"
    if command -v sha256sum >/dev/null 2>&1; then
        WISDOM_SHA256=$(sha256sum "$RESULT_DIR/fftw.wisdom" | awk '{print $1}')
    else
        WISDOM_SHA256=$(shasum -a 256 "$RESULT_DIR/fftw.wisdom" | awk '{print $1}')
    fi
    echo "fftw_wisdom_sha256=$WISDOM_SHA256" >> "$RESULT_DIR/run.status"
else
    echo "fftw_wisdom_sha256=not-generated" >> "$RESULT_DIR/run.status"
fi

first=true
for csv in "$RESULT_DIR"/concurrency_*.csv; do
    if [[ "$first" == true ]]; then
        cat "$csv" > "$COMBINED_CSV"
        first=false
    else
        tail -n +2 "$csv" >> "$COMBINED_CSV"
    fi
done

awk -F, '
    NR == 1 { next }
    {
        key = $1 SUBSEP $2 SUBSEP $3
        count[key]++
        raw[key] += ($6 == 0)
        sanitized[key] += ($7 == 0)
        plaintext[key] += ($8 == 0)
    }
    END {
        for (key in count) {
            split(key, fields, SUBSEP)
            print fields[1] "," fields[2] "," fields[3] "," count[key] "," \
                  raw[key] "," sanitized[key] "," plaintext[key]
        }
    }
' "$COMBINED_CSV" | sort -t, -k1,1 -k2,2n -k3,3n > "$OPERATION_CSV.tmp"

{
    echo "scope,threads,op,comparisons,raw_mismatches,sanitized_mismatches,plaintext_mismatches"
    cat "$OPERATION_CSV.tmp"
} > "$OPERATION_CSV"
rm -f "$OPERATION_CSV.tmp"

{
    echo "status=completed"
    echo "exit_code=0"
    echo "finished_at=$(date '+%F %T %z')"
    echo "summary=$COMBINED_CSV"
    echo "operation_summary=$OPERATION_CSV"
} >> "$RESULT_DIR/run.status"

echo "All raw ciphertext consensus checks passed."
echo "Results: $RESULT_DIR"
