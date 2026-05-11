#!/usr/bin/env python3
from __future__ import annotations

import os
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path

TARGET_OUTPUT = Path("/home/h2048/data/py/0416/myeloid_validation_optimized_cpu_rerun/adata_myeloid_refined_FINAL.h5ad")
OUTPUT_DIR = TARGET_OUTPUT.parent
SCVI_MODEL_DIR = OUTPUT_DIR / "scvi_model"
SCANVI_MODEL_DIR = OUTPUT_DIR / "scanvi_model"
PROCESS_PATTERN = "script/py/myeloid_scvi_scanvi_v2_4_cpu_rerun_20260416.py"
WRAPPER_CMD = [
    "/home/h2048/miniconda3/envs/scvi_env/bin/python",
    "/home/h2048/script/py/myeloid_scvi_scanvi_v2_4_cpu_rerun_20260416.py",
]
WORKDIR = Path("/home/h2048")
LOG_DIR = Path("/home/h2048/logs/20260417")
ORCH_LOG = LOG_DIR / "myeloid_scanvi_orchestrator_20260417.log"
RESTART_LOG = LOG_DIR / "myeloid_scanvi_restart_20260417.log"
CHECK_INTERVAL_SECONDS = 15
PROCESS_MISS_THRESHOLD = 2
MAX_RESTARTS = 3
CUTOVER_FLAG = OUTPUT_DIR / ".resumable_cutover_done"


def log(message: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    timestamped = f"[{datetime.now().strftime('%F %T')}] {message}"
    print(timestamped, flush=True)
    with ORCH_LOG.open("a", encoding="utf-8") as fh:
        fh.write(timestamped + "\n")


def find_matching_processes() -> list[tuple[int, str]]:
    proc = subprocess.run(["pgrep", "-af", PROCESS_PATTERN], capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not proc.stdout.strip():
        return []
    matches: list[tuple[int, str]] = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            pid_str, cmd = line.split(" ", 1)
            matches.append((int(pid_str), cmd))
        except ValueError:
            continue
    return matches


def restart_wrapper(reason: str) -> None:
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env.pop("PYTHONPATH", None)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with RESTART_LOG.open("a", encoding="utf-8") as log_fh:
        log_fh.write(f"\n[{datetime.now().strftime('%F %T')}] restarting wrapper: {reason}\n")
        log_fh.flush()
        subprocess.Popen(
            WRAPPER_CMD,
            cwd=str(WORKDIR),
            env=env,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    log(f"Wrapper restarted ({reason})")


def kill_processes(processes: list[tuple[int, str]], reason: str) -> None:
    target_pids: list[int] = []
    for pid, cmd in processes:
        try:
            os.kill(pid, signal.SIGTERM)
            log(f"Sent SIGTERM to PID {pid} ({reason}): {cmd}")
            target_pids.append(pid)
        except ProcessLookupError:
            log(f"PID {pid} already exited before SIGTERM ({reason})")

    if not target_pids:
        return

    time.sleep(10)
    for pid in target_pids:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            continue
        try:
            os.kill(pid, signal.SIGKILL)
            log(f"Escalated to SIGKILL for PID {pid} ({reason})")
        except ProcessLookupError:
            pass


def main() -> int:
    restart_count = 0
    consecutive_misses = 0
    last_proc_snapshot: tuple[tuple[int, str], ...] = ()

    log(f"Orchestrator started. Waiting for output: {TARGET_OUTPUT}")
    while True:
        if TARGET_OUTPUT.exists() and TARGET_OUTPUT.stat().st_size > 0:
            log(f"Target output detected; orchestrator exiting: {TARGET_OUTPUT}")
            return 0

        processes = find_matching_processes()
        proc_snapshot = tuple(processes)
        scvi_ready = SCVI_MODEL_DIR.exists()
        scanvi_ready = SCANVI_MODEL_DIR.exists()

        if scvi_ready and (not scanvi_ready) and (not CUTOVER_FLAG.exists()) and processes:
            log("Detected safe cutover point: scvi_model exists, scanvi_model absent, final output absent.")
            kill_processes(processes, reason="switch to resumable wrapper with scanvi early stopping")
            CUTOVER_FLAG.write_text("cutover_done\n", encoding="utf-8")
            time.sleep(5)
            restart_wrapper(reason="post-scvi cutover to resumable wrapper")
            consecutive_misses = 0
            last_proc_snapshot = ()
            time.sleep(CHECK_INTERVAL_SECONDS)
            continue

        if processes:
            consecutive_misses = 0
            if proc_snapshot != last_proc_snapshot:
                log("Active scanvi wrapper process(es):")
                for pid, cmd in processes:
                    log(f"  {pid} {cmd}")
                last_proc_snapshot = proc_snapshot
            time.sleep(CHECK_INTERVAL_SECONDS)
            continue

        consecutive_misses += 1
        log(
            f"No active scanvi wrapper process and no final output yet "
            f"(miss {consecutive_misses}/{PROCESS_MISS_THRESHOLD})."
        )
        if consecutive_misses < PROCESS_MISS_THRESHOLD:
            time.sleep(CHECK_INTERVAL_SECONDS)
            continue

        if restart_count >= MAX_RESTARTS:
            log("Maximum restart attempts reached without producing final output.")
            return 1

        restart_count += 1
        consecutive_misses = 0
        last_proc_snapshot = ()
        restart_wrapper(reason=f"failure recovery attempt {restart_count}/{MAX_RESTARTS}")
        time.sleep(CHECK_INTERVAL_SECONDS)


if __name__ == "__main__":
    raise SystemExit(main())
