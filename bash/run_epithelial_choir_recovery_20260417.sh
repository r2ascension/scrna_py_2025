#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/h2048"
R_SCRIPT="/home/h2048/script/R/epithelial_tissue_comparison_v1_3_2_choir_recovery_20260416.R"
R_BIN="/usr/bin/Rscript"
RUN_DATE="${RUN_DATE:-$(date '+%Y%m%d')}"
RUN_TS="${RUN_TS:-$(date '+%Y%m%d_%H%M%S')}"
LOG_DIR="/home/h2048/logs/${RUN_DATE}"
LOG_FILE="${LOG_DIR}/epithelial_choir_recovery_${RUN_TS}.log"
PID_FILE="${LOG_DIR}/epithelial_choir_recovery_${RUN_TS}.pid"
LATEST_PID_FILE="${LOG_DIR}/epithelial_choir_recovery_latest.pid"
LATEST_LOG_FILE="${LOG_DIR}/epithelial_choir_recovery_latest.log"
FORCE_START="${FORCE_START:-0}"
SCRIPT_PATTERN="${R_SCRIPT}"

mkdir -p "${LOG_DIR}"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

find_existing_pids() {
  pgrep -f "${SCRIPT_PATTERN}" || true
}

if [[ ! -f "${R_SCRIPT}" ]]; then
  log "[ERROR] Recovery R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

mapfile -t EXISTING_PIDS < <(find_existing_pids)
if (( ${#EXISTING_PIDS[@]} > 0 )) && [[ "${FORCE_START}" != "1" ]]; then
  log "[SKIP] Existing epithelial recovery process detected: ${EXISTING_PIDS[*]}"
  log "        Refusing to start a duplicate run."
  log "        If you really want to bypass this check, rerun with FORCE_START=1"
  exit 0
fi

cd "${ROOT_DIR}"

log "[START] Launching epithelial CHOIR recovery"
log "        R script : ${R_SCRIPT}"
log "        Log file : ${LOG_FILE}"
log "        PID file : ${PID_FILE}"

nohup /usr/bin/env -u LD_LIBRARY_PATH -u PYTHONPATH "${R_BIN}" "${R_SCRIPT}" >> "${LOG_FILE}" 2>&1 &
PID=$!
disown "${PID}" 2>/dev/null || true

printf '%s\n' "${PID}" > "${PID_FILE}"
printf '%s\n' "${PID}" > "${LATEST_PID_FILE}"
printf '%s\n' "${LOG_FILE}" > "${LATEST_LOG_FILE}"

log "[OK] Background process started: PID=${PID}"
log "[INFO] Tail log : tail -f ${LOG_FILE}"
log "[INFO] Check PID : ps -fp ${PID}"
log "[INFO] Stop run  : kill ${PID}"
