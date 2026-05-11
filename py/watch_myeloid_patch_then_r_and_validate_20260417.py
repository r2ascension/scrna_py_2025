#!/usr/bin/env python3
from __future__ import annotations

import ctypes
import os
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad

LIBC = ctypes.CDLL("libc.so.6", use_errno=True)
IN_CLOSE_WRITE = 0x00000008
IN_MOVED_TO = 0x00000080
WATCH_MASK = IN_CLOSE_WRITE | IN_MOVED_TO
EVENT_STRUCT_SIZE = 16

PATCH_OUTPUT = Path("/home/h2048/data/py/0416/adata_myeloid_L3refined_tissueaware_patched_v1.h5ad")
PATCH_DIR = PATCH_OUTPUT.parent
R_OUTPUT_DIR = Path("/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416")
FINAL_OUTPUTS = [
	PATCH_OUTPUT,
	R_OUTPUT_DIR / "myeloid_tissue_comparison_final.rds",
	R_OUTPUT_DIR / "myeloid_tissue_comparison_final.h5ad",
	R_OUTPUT_DIR / "REPORT.md",
]
R_CMD = [
	"/usr/bin/Rscript",
	"/home/h2048/script/R/myeloid_tissue_comparison_v1_2_3_20260416.R",
]
VALIDATE_CMD = [
	"/home/h2048/miniconda3/envs/scvi_env/bin/python",
	"/home/h2048/script/py/validate_myeloid_rerun_outputs_20260416.py",
]


def run_cmd(cmd: list[str]) -> None:
	env = os.environ.copy()
	env.pop("LD_LIBRARY_PATH", None)
	env.pop("PYTHONPATH", None)
	print(f"[patch-watcher] Running: {' '.join(cmd)}", flush=True)
	subprocess.run(cmd, check=True, env=env)


def print_outputs() -> None:
	for path in FINAL_OUTPUTS:
		status = "FOUND" if path.exists() else "MISSING"
		print(f"[patch-watcher] {status} {path}", flush=True)


def wait_for_stable_file(path: Path, checks: int = 3, delay_seconds: float = 2.0) -> None:
	stable_checks = 0
	previous_size: int | None = None
	while stable_checks < checks:
		if not path.exists():
			stable_checks = 0
			previous_size = None
			time.sleep(delay_seconds)
			continue

		current_size = path.stat().st_size
		if previous_size is not None and current_size == previous_size:
			stable_checks += 1
			print(
				f"[patch-watcher] Stable-size check {stable_checks}/{checks}: {path} ({current_size} bytes)",
				flush=True,
			)
		else:
			stable_checks = 0
			print(f"[patch-watcher] Waiting for size to stabilize: {path} ({current_size} bytes)", flush=True)
		previous_size = current_size
		time.sleep(delay_seconds)


def wait_until_readable(path: Path, delay_seconds: float = 2.0) -> None:
	while True:
		try:
			adata = ad.read_h5ad(path, backed="r")
			adata.file.close()
			print(f"[patch-watcher] Confirmed readable/unlocked: {path}", flush=True)
			return
		except BlockingIOError:
			print(f"[patch-watcher] File still locked, waiting: {path}", flush=True)
			time.sleep(delay_seconds)


def wait_for_patch_output() -> None:
	PATCH_DIR.mkdir(parents=True, exist_ok=True)

	if PATCH_OUTPUT.exists():
		print(f"[patch-watcher] Patch output already exists: {PATCH_OUTPUT}", flush=True)
		wait_for_stable_file(PATCH_OUTPUT)
		wait_until_readable(PATCH_OUTPUT)
		return

	print(f"[patch-watcher] Waiting for patch output: {PATCH_OUTPUT}", flush=True)

	fd = LIBC.inotify_init1(0)
	if fd < 0:
		err = ctypes.get_errno()
		raise OSError(err, os.strerror(err))
	try:
		wd = LIBC.inotify_add_watch(fd, os.fsencode(str(PATCH_DIR)), WATCH_MASK)
		if wd < 0:
			err = ctypes.get_errno()
			raise OSError(err, os.strerror(err))

		while not PATCH_OUTPUT.exists():
			data = os.read(fd, 4096)
			if not data:
				continue
			offset = 0
			while offset + EVENT_STRUCT_SIZE <= len(data):
				mask = ctypes.c_uint32.from_buffer_copy(data, offset + 4).value
				name_len = ctypes.c_uint32.from_buffer_copy(data, offset + 12).value
				offset += EVENT_STRUCT_SIZE
				name_bytes = data[offset:offset + name_len]
				offset += name_len
				name = name_bytes.split(b"\0", 1)[0].decode("utf-8", errors="ignore")
				if name:
					print(f"[patch-watcher] Event mask=0x{mask:x} name={name}", flush=True)
				if name == PATCH_OUTPUT.name and (mask & (IN_CLOSE_WRITE | IN_MOVED_TO)) and PATCH_OUTPUT.exists():
					print(f"[patch-watcher] Detected target file: {PATCH_OUTPUT}", flush=True)
					wait_for_stable_file(PATCH_OUTPUT)
					wait_until_readable(PATCH_OUTPUT)
					return
	finally:
		os.close(fd)


if __name__ == "__main__":
	try:
		wait_for_patch_output()
		run_cmd(R_CMD)
		run_cmd(VALIDATE_CMD)
		print_outputs()
		print("[patch-watcher] Downstream chain completed.", flush=True)
	except subprocess.CalledProcessError as exc:
		print(f"[patch-watcher] Command failed with exit code {exc.returncode}: {exc.cmd}", file=sys.stderr, flush=True)
		raise
	except Exception as exc:  # pragma: no cover
		print(f"[patch-watcher] Fatal error: {exc}", file=sys.stderr, flush=True)
		raise
