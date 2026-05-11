#!/usr/bin/env bash
set -u  # 保留未定义变量检查，不用 set -e 了

BASE_DIR="/home/h2048/data/source/1209/geo_downloads"
OUTPUT_ROOT="${BASE_DIR}/RAW_unpacked"

mkdir -p "${OUTPUT_ROOT}"

echo "Base dir     : ${BASE_DIR}"
echo "Output root  : ${OUTPUT_ROOT}"
echo "Searching for GSE*_RAW.tar ..."
echo

mapfile -t tar_files < <(find "${BASE_DIR}" -type f -name "GSE*_RAW.tar" | sort)

if [[ ${#tar_files[@]} -eq 0 ]]; then
    echo "没有找到任何 GSE*_RAW.tar 文件，检查路径是否正确。"
    exit 0
fi

for tar_file in "${tar_files[@]}"; do
    echo "Found: ${tar_file}"

    gse_dir="$(dirname "$(dirname "${tar_file}")")"
    gse_id="$(basename "${gse_dir}")"

    dest_dir="${OUTPUT_ROOT}/${gse_id}"
    mkdir -p "${dest_dir}"

    echo "  → GSE ID     : ${gse_id}"
    echo "  → Extract to : ${dest_dir}"

    if [[ -n "$(ls -A "${dest_dir}" 2>/dev/null || true)" ]]; then
        echo "  ⚠ 目标目录非空，认为已经解压过，跳过。"
        echo
        continue
    fi

    # 关键：允许 tar 报错，但脚本继续执行
    if tar -xf "${tar_file}" -C "${dest_dir}"; then
        echo "  ✓ 解压完成"
    else
        echo "  ✗ 解压失败，疑似文件损坏或下载不完整，跳过该 GSE。"
        # 可选：解压失败时清空目录
        rm -rf "${dest_dir}"
    fi
    echo
done

echo "全部完成。"
