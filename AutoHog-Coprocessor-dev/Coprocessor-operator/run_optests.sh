#!/bin/bash
# bash ./run_optests.sh
# ========================================================

# ========================================================

EXE_PATH="../build/Coprocessor-operator/main_testop" 

# ========================================================

# ========================================================
BIT_WIDTH=8
MODE="triv"
EPOCHS=100
DETAIL=0


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

LOG_DIR="${SCRIPT_DIR}/fhe_operations_logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/fhe_test_${TIMESTAMP}.log"

# ========================================================

# ========================================================


if [ ! -d "$LOG_DIR" ]; then
    echo "📂 创建日志文件夹: $LOG_DIR"
    mkdir -p "$LOG_DIR"
fi


if [ ! -f "$EXE_PATH" ]; then
    echo "❌ 错误: 找不到可执行文件 $EXE_PATH"
    exit 1
fi


OPERATIONS=(
    "FheAdd" "FheSub" 
    "FheBitAnd" "FheBitOr" "FheBitXor" "FheBitNot" 
    "FheEq" "FheNe" "FheGe" "FheGt" "FheLe" "FheLt" 
    "FheIfThenElse" "TrivialEncrypt"
)

# OPERATIONS=("FheIfThenElse" )



echo "🚀 开始自动化同态算子测试..."
echo "📂 日志目录: $LOG_DIR"
echo "📄 当前日志: $(basename "$LOG_FILE")"
echo "--------------------------------------------------------"


{
    echo "FHE OPERATOR AUTOMATION TEST REPORT"
    echo "Date: $(date)"
    echo "Config: Bits=${BIT_WIDTH}, Mode=${MODE}, Epochs=${EPOCHS}"
    echo "========================================================"
} > "$LOG_FILE"


for op in "${OPERATIONS[@]}"
do
    echo -n "正在执行: [$op] ... "
    
    echo -e "\n\n>>> OPERATION: $op" >> "$LOG_FILE"
    echo "--------------------------------------------------------" >> "$LOG_FILE"
    

    "$EXE_PATH" "$op" "$BIT_WIDTH" "$MODE" "$EPOCHS" "$DETAIL" >> "$LOG_FILE" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✅ 完成"
    else
        echo "❌ 出错 (请查阅日志)"
    fi
done

echo "--------------------------------------------------------"
echo "🎉 所有测试任务已跑完！"
echo "📊 结果汇总请查看: $LOG_FILE"