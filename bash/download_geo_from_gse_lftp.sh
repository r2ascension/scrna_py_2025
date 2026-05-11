#!/usr/bin/env bash
# 批量从 GEO FTP 下载 GSE 的 suppl 和 matrix 文件（lftp + 每文件进度）
# 用法:
#   ./download_geo_from_gse_lftp_pget.sh gse_list.txt /path/to/output
#
# 依赖:
#   - lftp (你的是 4.9.2，OK)

set -u  # 使用未定义变量时报错

GSE_LIST_FILE="${1:-gse_list.txt}"
OUTDIR="${2:-./geo_downloads}"

FTP_HOST="ftp.ncbi.nlm.nih.gov"
FTP_BASE="/geo/series"

if [ ! -f "$GSE_LIST_FILE" ]; then
    echo "找不到 GSE 列表文件: $GSE_LIST_FILE"
    exit 1
fi

mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# 工具函数：用 lftp+pget 下载某个远程目录中的所有文件
#   参数1：远程目录（例如 /geo/series/GSE155nnn/GSE155113/suppl）
#   参数2：本地目录（例如 ./geo_downloads/GSE155113/suppl）
# ------------------------------------------------------------
download_dir_with_pget() {
    local remote_dir="$1"
    local local_dir="$2"

    mkdir -p "$local_dir"

    echo "  - 列出远程目录: $remote_dir"

    # 列出远程目录下的文件（可能是文件名，也可能是带路径，统一用 basename 处理）
    local files
    if ! files=$(lftp -u anonymous, "$FTP_HOST" -e "cls -1 \"$remote_dir\"; bye" 2>/dev/null); then
        echo "    (远程目录不存在或无法访问，跳过)"
        return 0
    fi

    if [ -z "$files" ]; then
        echo "    (目录为空，跳过)"
        return 0
    fi

    echo "$files" | while IFS= read -r fname; do
        [ -z "$fname" ] && continue

        # 防止 cls 返回的是完整路径，这里统一只用文件名部分
        local base
        base="$(basename "$fname")"

        echo "    ▶ 下载文件: $base"

        lftp -u anonymous, "$FTP_HOST" <<EOF
set cmd:show-status yes
set cmd:status-interval 1
set xfer:rate-period 5
set xfer:eta-period 5
set xfer:eta-terse no

set net:timeout 20
set net:max-retries 20
set net:reconnect-interval-base 5
set net:reconnect-interval-max 60
set ftp:passive-mode yes

# 每个文件单独 pget，支持多线程 + 断点续传
pget -n 4 -c "$remote_dir/$base" -o "$local_dir/$base"

bye
EOF

    done
}

echo ">>> 使用列表文件: $GSE_LIST_FILE"
echo ">>> 输出目录: $OUTDIR"
echo

# 主循环：逐个 GSE 处理
while read -r gse_raw; do
    # 去掉空白
    gse="$(echo "$gse_raw" | tr -d '[:space:]')"
    [ -z "$gse" ] && continue

    # 统一大写
    gse="$(echo "$gse" | tr '[:lower:]' '[:upper:]')"

    # 必须以 GSE 开头
    if [[ ! "$gse" =~ ^GSE[0-9]+$ ]]; then
        echo "!!! 跳过无效 GSE ID: $gse_raw"
        continue
    fi

    num="${gse#GSE}"  # 去掉前缀 GSE，剩下纯数字

    # GEO FTP 目录规则：
    #   GSE999      -> /geo/series/GSEnnn/GSE999/...
    #   GSE1000     -> /geo/series/GSE1nnn/GSE1000/...
    #   GSE50499    -> /geo/series/GSE50nnn/GSE50499/...
    if [ "${#num}" -gt 3 ]; then
        series_prefix="${num::-3}"    # 去掉最后 3 位
        series_dir="GSE${series_prefix}nnn"
    else
        series_dir="GSEnnn"
    fi

    echo "================================================================"
    echo "处理: $gse (series 目录: $series_dir)"
    echo "================================================================"

    # 本地目录
    gse_out="$OUTDIR/$gse"
    local_suppl="$gse_out/suppl"
    local_matrix="$gse_out/matrix"
    mkdir -p "$local_suppl" "$local_matrix"

    # 远程目录
    remote_suppl="$FTP_BASE/$series_dir/$gse/suppl"
    remote_matrix="$FTP_BASE/$series_dir/$gse/matrix"

    echo "[1/2] suppl 目录:  $remote_suppl"
    download_dir_with_pget "$remote_suppl" "$local_suppl"

    echo "[2/2] matrix 目录: $remote_matrix"
    download_dir_with_pget "$remote_matrix" "$local_matrix"

    echo
done < "$GSE_LIST_FILE"

echo ">>> 全部 GSE 处理完成。文件保存在: $OUTDIR"
