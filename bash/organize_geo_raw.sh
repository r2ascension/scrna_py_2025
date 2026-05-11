#!/usr/bin/env bash
# 批量解压 GEO GSE*_RAW.tar，并按 GSM 归类到子目录
# 用法:
#   ./organize_geo_raw.sh /path/to/raw_tar_dir /data/GEO
#
# 参数:
#   $1 = 包含 GSE*_RAW.tar 的目录（Linux 路径）
#   $2 = DATA_ROOT，对应 R 脚本中的 DATA_ROOT（默认 /data/GEO）

set -euo pipefail

if [[ "${1:-}" == "" || "${2:-}" == "" ]]; then
    echo "用法: $0 RAW_TAR_DIR DATA_ROOT"
    echo "例如:"
    echo "  $0 /mnt/e/idm /data/GEO"
    exit 1
fi

RAW_TAR_DIR="$1"
DATA_ROOT="$2"

if [[ ! -d "$RAW_TAR_DIR" ]]; then
    echo "错误: RAW_TAR_DIR 不存在: $RAW_TAR_DIR"
    exit 1
fi

mkdir -p "$DATA_ROOT"

# 确保通配符匹配不到文件时不会原样展开
shopt -s nullglob

echo "============================================================"
echo "GEO RAW.tar 解压 & GSM 分组脚本"
echo "RAW_TAR_DIR: $RAW_TAR_DIR"
echo "DATA_ROOT  : $DATA_ROOT"
echo "============================================================"
echo

# ------------------------------------------------------------------
# Step 1: 解压每个 GSE*_RAW.tar 到对应的 GSE 目录
# ------------------------------------------------------------------

tars=( "$RAW_TAR_DIR"/GSE*_RAW.tar "$RAW_TAR_DIR"/GSE*_RAW.tar.gz )
if (( ${#tars[@]} == 0 )); then
    echo "在 $RAW_TAR_DIR 下未找到 GSE*_RAW.tar 或 GSE*_RAW.tar.gz 文件"
    exit 1
fi

for tarfile in "${tars[@]}"; do
    bname="$(basename "$tarfile")"

    # 从文件名中提取 GSE ID（例如 GSE202100_RAW.tar → GSE202100）
    gse_id="${bname%%_*}"  # 取第一个下划线前面的部分
    if [[ ! "$gse_id" =~ ^GSE[0-9]+$ ]]; then
        echo "⚠️  跳过无法识别 GSE ID 的文件: $bname"
        continue
    fi

    gse_dir="$DATA_ROOT/$gse_id"
    mkdir -p "$gse_dir"

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "处理 $bname"
    echo "  GSE ID    : $gse_id"
    echo "  解压目录  : $gse_dir"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # 解压到对应 GSE 目录
    tar -xf "$tarfile" -C "$gse_dir"

    echo "  ✓ 解压完成: $tarfile → $gse_dir"
    echo
done

# ------------------------------------------------------------------
# Step 2: 在每个 GSE 目录中按 GSM 归类文件
#   规则:
#   - 所有路径中包含 'GSM[0-9]+' 的文件
#   - 将其移动到 $DATA_ROOT/GSExxxx/GSMyyyyyyy/ 下
#   - 已经在 GSM 目录中的文件会自动跳过（做到幂等）
# ------------------------------------------------------------------

echo
echo "============================================================"
echo "按 GSM 归类文件..."
echo "============================================================"

# 遍历 DATA_ROOT 下所有 GSE 目录
for gse_dir in "$DATA_ROOT"/GSE[0-9]*; do
    [[ -d "$gse_dir" ]] || continue
    gse_id="$(basename "$gse_dir")"
    echo
    echo "▶ 处理 GSE 目录: $gse_id ($gse_dir)"

    # 查找该 GSE 目录下的所有文件（包含子目录），逐个处理
    # 使用 -print0 以安全处理特殊字符
    while IFS= read -r -d '' file; do
        # file 的完整路径，例如 /data/GEO/GSE202100/GSM4695772_xxx.txt.gz
        dir_name="$(basename "$(dirname "$file")")"

        # 如果已经在 GSM 目录内（父目录名是 GSMxxxxxx），跳过
        if [[ "$dir_name" =~ ^GSM[0-9]+$ ]]; then
            continue
        fi

        fname="$(basename "$file")"

        # 在文件名中匹配 GSM ID
        if [[ "$fname" =~ (GSM[0-9]+) ]]; then
            gsm_id="${BASH_REMATCH[1]}"
        else
            echo "  ⚠️  未找到 GSM ID，跳过: $fname"
            continue
        fi

        gsm_dir="$gse_dir/$gsm_id"
        mkdir -p "$gsm_dir"

        echo "  移动文件: $fname → $gsm_id/"
        mv "$file" "$gsm_dir/"

    done < <(find "$gse_dir" -type f -print0)

    echo "  ✓ $gse_id 处理完成"
done

echo
echo "============================================================"
echo "✅ 所有 GSE_RAW.tar 解压完成，并按 GSM 归类"
echo "ROOT 目录结构示例:"
echo "  $DATA_ROOT/GSE155113/GSMxxxxxx/*"
echo "  $DATA_ROOT/GSE156285/GSMyyyyyy/*"
echo "============================================================"
