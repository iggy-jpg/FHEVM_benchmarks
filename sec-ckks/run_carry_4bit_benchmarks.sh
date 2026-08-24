#!/usr/bin/env bash

# Run the 4-bit multiplication benchmarks and keep one log/status pair per
# precision. Set CARRY_JOBS to control how many precisions run concurrently.
set -u

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BIN="${SCRIPT_DIR}/build/bin/examples/pke/carry_test"
OUT_DIR="${SCRIPT_DIR}/build/testdata/carry"
CARRY_JOBS=${CARRY_JOBS:-2}

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    echo "Build it with: cmake --build ${SCRIPT_DIR}/build --target carry_test -j2" >&2
    exit 1
fi

if ! [[ "${CARRY_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "CARRY_JOBS must be a positive integer" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

if (($# == 0)); then
    set -- 8 16 32 64 128 256
fi

run_one() {
    local bits=$1
    local chunks
    local pbits=8
    local p
    local log
    local status
    local rc

    if ! [[ "${bits}" =~ ^[1-9][0-9]*$ ]] || ((bits % 4 != 0)); then
        echo "Skipping invalid precision: ${bits}" >&2
        return 0
    fi
    chunks=$((bits / 4))
    if ((chunks > 64)); then
        echo "Skipping ${bits}-bit: segmented carry currently supports at most 256-bit inputs" >&2
        return 0
    fi

    p=$((1 << pbits))
    log="${OUT_DIR}/mul_${bits}bit_4bit_p${p}.log"
    status="${OUT_DIR}/mul_${bits}bit_4bit_p${p}.status"

    {
        echo "Command: ${BIN} ${bits} 4 mul"
        echo "Bootstrap bits: ${pbits} (p=${p})"
        echo "Started: $(date --iso-8601=seconds)"
        /usr/bin/time -v "${BIN}" "${bits}" 4 mul
        rc=$?
        echo "Finished: $(date --iso-8601=seconds)"
        echo "Exit status: ${rc}"
    } >"${log}" 2>&1

    printf '%s\n' "${rc}" >"${status}"
    echo "${bits}-bit: status=${rc}, log=${log}"
}

pending=()
for bits in "$@"; do
    run_one "${bits}" &
    pending+=("$!")
    if ((${#pending[@]} >= CARRY_JOBS)); then
        wait "${pending[0]}" || true
        pending=("${pending[@]:1}")
    fi
done

for pid in "${pending[@]}"; do
    wait "${pid}" || true
done
