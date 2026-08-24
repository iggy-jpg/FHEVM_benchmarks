#!/usr/bin/env bash
set -euo pipefail

# Compare test_with_zk under:
#   1. limited resources: 4 CPU cores, 8 GB RAM by default
#   2. normal resources: no script-imposed CPU/RAM limits
#
# Examples:
#   ./run_test_with_zk.sh
#   NUM_BITS=256 ./run_test_with_zk.sh
#   CRS_FILE=../zk_crs.bin ./run_test_with_zk.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build"
BINARY="${BUILD_DIR}/Coprocessor-operator/test_with_zk"
LOG_ROOT="${SCRIPT_DIR}/fhe_operations_logs"
RUN_ID="$(date +"%Y%m%d_%H%M%S")"
COMPARE_DIR="${LOG_ROOT}/test_with_zk_compare_${RUN_ID}"

CPU_CORES="${CPU_CORES:-4}"
MEM_LIMIT_GB="${MEM_LIMIT_GB:-8}"
NUM_BITS="${NUM_BITS:-128}"

MEM_LIMIT_KB=$((MEM_LIMIT_GB * 1024 * 1024))

echo "=========================================================="
echo "Comparing test_with_zk"
echo "Project root: ${PROJECT_ROOT}"
echo "Limited CPU:  ${CPU_CORES} cores"
echo "Limited RAM:  ${MEM_LIMIT_GB} GB"
echo "Payload bits: ${NUM_BITS}"
echo "=========================================================="

if [ ! -x "${BINARY}" ]; then
    echo "test_with_zk binary not found; building target..."
    cmake -S "${PROJECT_ROOT}" -B "${BUILD_DIR}"
    cmake --build "${BUILD_DIR}" --target test_with_zk --parallel "${CPU_CORES}"
fi

if [ -n "${CRS_FILE:-}" ]; then
    CRS_PATH="${CRS_FILE}"
elif [ -f "${PROJECT_ROOT}/crs_lvl2.bin" ]; then
    CRS_PATH="${PROJECT_ROOT}/crs_lvl2.bin"
elif [ -f "${SCRIPT_DIR}/crs_lvl2.bin" ]; then
    CRS_PATH="${SCRIPT_DIR}/crs_lvl2.bin"
else
    CRS_PATH="${PROJECT_ROOT}/crs_lvl2.bin"
    GEN_CRS="${BUILD_DIR}/Coprocessor-operator/gen_crs"
    if [ ! -x "${GEN_CRS}" ]; then
        echo "gen_crs binary not found; building target..."
        cmake -S "${PROJECT_ROOT}" -B "${BUILD_DIR}"
        cmake --build "${BUILD_DIR}" --target gen_crs --parallel "${CPU_CORES}"
    fi
    echo "CRS file not found; generating ${CRS_PATH} ..."
    "${GEN_CRS}" "${CRS_PATH}"
fi

CRS_DIGEST="$(shasum -a 256 "${CRS_PATH}" | awk '{print $1}')"

echo "CRS file:     ${CRS_PATH}"
echo "CRS sha256:   ${CRS_DIGEST}"
echo "Executable:   ${BINARY}"
echo "Log dir:      ${COMPARE_DIR}"
echo "----------------------------------------------------------"

CMD=("${BINARY}" "${NUM_BITS}" "--crs" "${CRS_PATH}" "--crs-sha256" "${CRS_DIGEST}")

mkdir -p "${COMPARE_DIR}"

run_limited() {
    local log_file="$1"

    (
        export OMP_NUM_THREADS="${CPU_CORES}"
        export OPENBLAS_NUM_THREADS="${CPU_CORES}"
        export MKL_NUM_THREADS="${CPU_CORES}"
        export VECLIB_MAXIMUM_THREADS="${CPU_CORES}"
        export NUMEXPR_NUM_THREADS="${CPU_CORES}"

        if ! ulimit -v "${MEM_LIMIT_KB}" 2>/dev/null; then
            echo "WARN: failed to set virtual memory limit with ulimit -v."
        fi

        if ! ulimit -d "${MEM_LIMIT_KB}" 2>/dev/null; then
            echo "WARN: failed to set data segment memory limit with ulimit -d."
        fi

        if command -v taskset >/dev/null 2>&1; then
            LAST_CORE=$((CPU_CORES - 1))
            exec taskset -c "0-${LAST_CORE}" "${CMD[@]}"
        elif command -v cpulimit >/dev/null 2>&1; then
            exec cpulimit -l "$((CPU_CORES * 100))" -- "${CMD[@]}"
        else
            echo "WARN: taskset/cpulimit not found; CPU is limited through thread-count env vars only."
            exec "${CMD[@]}"
        fi
    ) >"${log_file}" 2>&1
}

run_normal() {
    local log_file="$1"

    (
        exec "${CMD[@]}"
    ) >"${log_file}" 2>&1
}

run_case() {
    local name="$1"
    local log_file="$2"

    echo ">>> Running ${name} case..."
    local start_ts
    start_ts="$(date +%s)"

    set +e
    if [ "${name}" = "limited" ]; then
        run_limited "${log_file}"
    else
        run_normal "${log_file}"
    fi
    local status=$?
    set -e

    local end_ts
    end_ts="$(date +%s)"
    local duration=$((end_ts - start_ts))

    printf "%s,%s,%s\n" "${name}" "${status}" "${duration}" >> "${COMPARE_DIR}/summary.csv"

    if [ "${status}" -eq 0 ]; then
        echo "    ${name}: PASS (${duration}s)"
    else
        echo "    ${name}: FAIL status=${status} (${duration}s)"
        FAILED_CASES=$((FAILED_CASES + 1))
    fi
    echo "    log: ${log_file}"
}

echo "case,status,duration_seconds" > "${COMPARE_DIR}/summary.csv"
FAILED_CASES=0

run_case "limited" "${COMPARE_DIR}/limited.log"
run_case "normal" "${COMPARE_DIR}/normal.log"

echo "----------------------------------------------------------"
echo "Comparison finished."
echo "Summary: ${COMPARE_DIR}/summary.csv"

if [ "${FAILED_CASES}" -ne 0 ]; then
    exit 1
fi
