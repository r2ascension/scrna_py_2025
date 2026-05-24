#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

TARGET_OUTPUT = Path("/home/h2048/data/py/0416/myeloid_validation_optimized_cpu_rerun/adata_myeloid_refined_FINAL.h5ad")
PROCESS_PATTERN = "script/py/myeloid_scvi_scanvi_v2_4_cpu_rerun_20260416.py"
WRAPPER_CMD = [
    "/home/h2048/miniconda3/envs/scvi_env/bin/python",
    "/home/h2048/script/py/myeloid_scvi_scanvi_v2_4_cpu_rerun_20260416.py",
]
WORKDIR = Path("/home/h2048")
LOG_DIR = Path("/home/h2048/logs/20260417")
SUPERVISOR_LOG = LOG_DIR / "myeloid_scanvi_supervisor_20260417.log"
RESTART_LOG = LOG_DIR / "myeloid_scanvi_restart_20260417.log"
CHECK_INTERVAL_SECONDS = 120
PROCESS_MISS_THRESHOLD = 2
MAX_RESTARTS = 3


def log(message: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    timestamped = f"[{datetime.now().strftime('%F %T')}] {message}"
    print(timestamped, flush=True)
    with SUPERVISOR_LOG.open("a", encoding="utf-8") as fh:
        fh.write(timestamped + "\n")


def find_matching_processes() -> list[str]:
    proc = subprocess.run(
        ["pgrep", "-af", PROCESS_PATTERN],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        return []
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def restart_wrapper() -> None:
    env = os.environ.copy()
    env.pop("LD_LIBRARY_PATH", None)
    env.pop("PYTHONPATH", None)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with RESTART_LOG.open("a", encoding="utf-8") as log_fh:
        log_fh.write(f"\n[{datetime.now().strftime('%F %T')}] restarting wrapper\n")
        log_fh.flush()
        subprocess.Popen(
            WRAPPER_CMD,
            cwd=str(WORKDIR),
            env=env,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )


def main() -> int:
    restart_count = 0
    consecutive_misses = 0
    last_proc_snapshot: tuple[str, ...] = ()

    log(f"Supervisor started. Waiting for output: {TARGET_OUTPUT}")
    while True:
        if TARGET_OUTPUT.exists() and TARGET_OUTPUT.stat().st_size > 0:
            log(f"Target output detected; supervisor exiting: {TARGET_OUTPUT}")
            return 0

        processes = find_matching_processes()
        proc_snapshot = tuple(processes)
        if processes:
            consecutive_misses = 0
            if proc_snapshot != last_proc_snapshot:
                log("Active scanvi wrapper process(es):")
                for line in processes:
                    log(f"  {line}")
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
        log(f"Restarting CPU-safe myeloid scanvi wrapper (attempt {restart_count}/{MAX_RESTARTS})...")
        restart_wrapper()
        time.sleep(CHECK_INTERVAL_SECONDS)


if __name__ == "__main__":
    raise SystemExit(main())
