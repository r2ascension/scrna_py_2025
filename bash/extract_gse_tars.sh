#!/usr/bin/env bash
# 解压指定 GSE 目录下各 GSM 子目录中的 *.tar.gz 到其所在目录，并删除源压缩包
# 用法:
#   ./extract_gse_tars.sh /home/h2048/data/source/1210/GSE235711

set -euo pipefail

GSE_DIR="${1:-}"

if [[ -z "$GSE_DIR" ]]; then
    echo "用法: $0 /path/to/GSE_dir"
    exit 1
fi

if [[ ! -d "$GSE_DIR" ]]; then
    echo "错误: 目录不存在: $GSE_DIR"
    exit 1
fi

echo ">>> 处理目录: $GSE_DIR"
shopt -s nullglob

# 遍历 GSE 下所有 GSM 子目录里的 tar.gz
for tarfile in "$GSE_DIR"/GSM*/*.tar.gz; do
    [ -f "$tarfile" ] || continue
    gsm_dir="$(dirname "$tarfile")"

    echo "------------------------------------------------------------"
    echo "GSM 目录: $gsm_dir"
    echo "压缩包:   $(basename "$tarfile")"

    # 解压到当前 GSM 目录
    tar -xzf "$tarfile" -C "$gsm_dir"

    # 仅在解压成功后删除
    rm -v "$tarfile"
done

echo ">>> 完成"
