#!/usr/bin/env bash
# 批量从 GEO FTP 下载 GSE 的 suppl 和 matrix 文件（lftp + 每文件校验）
# 用法:
#   ./download_geo_from_gse_lftp_pget.sh gse_list.txt /path/to/output
#
# 依赖:
#   - lftp

set -u  # 使用未定义变量时报错

GSE_LIST_FILE="${1:-gse_list.txt}"
OUTDIR="${2:-./geo_downloads}"

FTP_HOST="ftp.ncbi.nlm.nih.gov"
FTP_BASE="/geo/series"

# 下载失败 / 校验失败的日志
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${OUTDIR}" "${LOG_DIR}"
FAIL_LOG="${LOG_DIR}/download_failures_$(date +%Y%m%d_%H%M%S).log"
touch "${FAIL_LOG}"

if [ ! -f "$GSE_LIST_FILE" ]; then
    echo "找不到 GSE 列表文件: $GSE_LIST_FILE"
    exit 1
fi

# ------------------------------------------------------------
# 工具函数：对单个文件做完整性检查
#   返回 0 = 通过；1 = 失败
# ------------------------------------------------------------
verify_file() {
    local file="$1"

    if [[ ! -f "$file" ]]; then
        echo "      [check] 文件不存在: $file"
        return 1
    fi

    if [[ ! -s "$file" ]]; then
        echo "      [check] 文件大小为 0: $file"
        return 1
    fi

    case "$file" in
        *.tar|*.tar.gz)
            # 检查 tar 结构是否正常（截断会报错）
            if tar -tf "$file" >/dev/null 2>&1; then
                echo "      [check] tar 结构检查通过"
                return 0
            else
                echo "      [check] tar 结构检查失败 (可能截断/损坏)"
                return 1
            fi
            ;;
        *.gz)
            # 检查 gzip 完整性——GEO 官方也推荐用 gunzip -t 检查 gz 文件:contentReference[oaicite:1]{index=1}
            if gzip -t "$file" >/dev/null 2>&1; then
                echo "      [check] gzip 完整性检查通过"
                return 0
            else
                echo "      [check] gzip 完整性检查失败 (压缩文件损坏)"
                return 1
            fi
            ;;
        *)
            # 其他类型目前只能做“非空”检查
            echo "      [check] 非压缩文件，跳过深度校验，仅检测非 0 字节"
            return 0
            ;;
    esac
}

# ------------------------------------------------------------
# 工具函数：用 lftp+pget 下载“目录内的所有文件”（逐文件 + 校验 + 重试）
#   参数1：远程目录（例如 /geo/series/GSE155nnn/GSE155113/suppl）
#   参数2：本地目录（例如 ./geo_downloads/GSE155113/suppl）
# ------------------------------------------------------------
download_dir_with_pget() {
    local remote_dir="$1"
    local local_dir="$2"

    mkdir -p "$local_dir"

    echo "  - 列出远程目录: $remote_dir"

    # 列出远程目录下的文件名（可能是文件名，也可能是带路径，统一用 basename）
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

        local base
        base="$(basename "$fname")"
        local remote_file="${remote_dir}/${base}"
        local local_file="${local_dir}/${base}"

        # 已存在且非空可以选择跳过（防止反复重下）
        if [[ -s "$local_file" ]]; then
            echo "    ▶ 本地已存在非空文件，先做校验: $base"
            if verify_file "$local_file"; then
                echo "    ✓ 已有文件通过校验，跳过下载"
                continue
            else
                echo "    ⚠ 已有文件校验失败，将重新下载"
                rm -f "$local_file"
            fi
        fi

        # 最多重试次数
        local max_attempts=3
        local attempt
        local success=0

        for ((attempt=1; attempt<=max_attempts; attempt++)); do
            echo "    ▶ 下载文件 (尝试 ${attempt}/${max_attempts}): $base"

            lftp -u anonymous, "$FTP_HOST" <<EOF
set cmd:show-status yes
set cmd:status-interval 5
set xfer:rate-period 5
set xfer:eta-period 5
set xfer:eta-terse no

set net:timeout 60
set net:max-retries 20
set net:reconnect-interval-base 5
set net:reconnect-interval-max 60
set ftp:passive-mode yes

# 使用 pget 支持多线程 + 断点续传
pget -n 4 -c "$remote_file" -o "$local_file"

bye
EOF

            echo "      - 下载完成，开始完整性检查..."
            if verify_file "$local_file"; then
                echo "    ✓ 文件下载并通过校验: $base"
                success=1
                break
            else
                echo "    ✗ 校验失败，将删除并重试: $base"
                rm -f "$local_file"
            fi
        done

        if [[ $success -ne 1 ]]; then
            echo "    !!! 多次重试后仍失败，记录日志并跳过: $base"
            echo "$(date '+%F %T') | ${remote_file} | verify_failed" >> "${FAIL_LOG}"
        fi

        echo
    done
}

echo ">>> 使用列表文件: $GSE_LIST_FILE"
echo ">>> 输出目录: $OUTDIR"
echo ">>> 失败日志: $FAIL_LOG"
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
echo ">>> 若有失败文件，请查看: $FAIL_LOG"
