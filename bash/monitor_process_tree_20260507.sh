#!/usr/bin/env bash
set -u

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <root_pid> <output_tsv> [interval_seconds] [label]" >&2
  exit 2
fi

ROOT_PID="$1"
OUT_TSV="$2"
INTERVAL="${3:-30}"
LABEL="${4:-process}"

if ! [[ "$ROOT_PID" =~ ^[0-9]+$ ]]; then
  echo "root_pid must be numeric: $ROOT_PID" >&2
  exit 2
fi
if ! [[ "$INTERVAL" =~ ^[0-9]+$ ]] || [[ "$INTERVAL" -lt 1 ]]; then
  INTERVAL=30
fi

mkdir -p "$(dirname "$OUT_TSV")"
if [[ ! -s "$OUT_TSV" ]]; then
  printf 'timestamp\tlabel\troot_pid\tphase\tpids\tn_proc\trss_kb_sum\tvsz_kb_sum\tpcpu_sum\tpmem_sum\tmem_total_kb\tmem_available_kb\ttop_processes\n' > "$OUT_TSV"
fi

collect_pids() {
  /usr/bin/python3 - "$ROOT_PID" <<'PY'
import subprocess
import sys
from collections import defaultdict, deque

root = int(sys.argv[1])
try:
    raw = subprocess.check_output(["ps", "-eo", "pid=,ppid="], text=True)
except Exception:
    print(root)
    raise SystemExit
children = defaultdict(list)
seen_pids = set()
for line in raw.splitlines():
    parts = line.split()
    if len(parts) < 2:
        continue
    try:
        pid = int(parts[0]); ppid = int(parts[1])
    except ValueError:
        continue
    seen_pids.add(pid)
    children[ppid].append(pid)

if root not in seen_pids:
    print(root)
    raise SystemExit

seen = []
visited = set()
queue = deque([root])
while queue:
    pid = queue.popleft()
    if pid in visited:
        continue
    visited.add(pid)
    seen.append(pid)
    for child in children.get(pid, []):
        queue.append(child)
print(",".join(str(x) for x in seen))
PY
}

sample_once() {
  local phase="$1"
  local now pids ps_rows n_proc rss_sum vsz_sum pcpu_sum pmem_sum mem_total mem_avail top_processes
  now="$(date '+%Y-%m-%d %H:%M:%S')"
  pids="$(collect_pids 2>/dev/null || true)"
  if [[ -z "$pids" ]]; then
    pids="$ROOT_PID"
  fi

  ps_rows="$(ps -o pid=,ppid=,stat=,etime=,pcpu=,pmem=,rss=,vsz=,comm= -p "$pids" 2>/dev/null || true)"
  n_proc="$(printf '%s\n' "$ps_rows" | awk 'NF {n++} END {print n+0}')"
  rss_sum="$(printf '%s\n' "$ps_rows" | awk 'NF {s+=$7} END {printf "%.0f", s+0}')"
  vsz_sum="$(printf '%s\n' "$ps_rows" | awk 'NF {s+=$8} END {printf "%.0f", s+0}')"
  pcpu_sum="$(printf '%s\n' "$ps_rows" | awk 'NF {s+=$5} END {printf "%.2f", s+0}')"
  pmem_sum="$(printf '%s\n' "$ps_rows" | awk 'NF {s+=$6} END {printf "%.2f", s+0}')"
  mem_total="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null || echo NA)"
  mem_avail="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo 2>/dev/null || echo NA)"
  top_processes="$(printf '%s\n' "$ps_rows" | awk 'NF {printf "%s:%sKB:%s;", $1, $7, $9}' | sed 's/[[:space:]\t]/_/g')"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$now" "$LABEL" "$ROOT_PID" "$phase" "$pids" "$n_proc" "$rss_sum" "$vsz_sum" \
    "$pcpu_sum" "$pmem_sum" "$mem_total" "$mem_avail" "$top_processes" >> "$OUT_TSV"
}

sample_once "start"
while kill -0 "$ROOT_PID" 2>/dev/null; do
  sample_once "running"
  sleep "$INTERVAL"
done
sample_once "finished"
