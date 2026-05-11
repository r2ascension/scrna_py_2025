#!/usr/bin/env python3
"""CPU-safe wrapper for the validated myeloid scVI/scANVI production script.

Purpose:
- keep `myeloid_scvi_scanvi_v2_4_20260322.py` immutable
- force CPU execution because the current host NVIDIA driver is too old for the
  installed PyTorch CUDA build
- rerun into a versioned 2026-04-16 output directory instead of overwriting the
  historical 2026-03-22 production outputs
"""

import os
import re
from pathlib import Path

BASE_SCRIPT = Path("/home/h2048/script/py/myeloid_scvi_scanvi_v2_4_20260322.py")
OLD_OUTPUT_LINE = 'OUTPUT_DIR = Path("/home/h2048/data/py/0322/myeloid_validation_optimized")'
NEW_OUTPUT_LINE = 'OUTPUT_DIR = Path("/home/h2048/data/py/0416/myeloid_validation_optimized_cpu_rerun")'
SCVI_BLOCK_PATTERN = re.compile(
  r'print\("Training scVI\.\.\."\)\n'
  r'vae\.train\(max_epochs=MAX_EPOCHS_SCVI, batch_size=BATCH_SIZE, early_stopping=True\)\n'
  r'print\("  Training complete"\)\n\n'
  r'vae\.save\(OUTPUT_DIR / "scvi_model", overwrite=True\)\n'
  r'print\("  Model saved"\)',
  re.MULTILINE,
)
SCANVI_BLOCK_PATTERN = re.compile(
  r'print\("Training scANVI\.\.\."\)\n'
  r'lvae = scvi\.model\.SCANVI\.from_scvi_model\(\n'
  r'\s+vae, adata=adata_scvi, labels_key=CELLTYPE_L3_REFINED_COL, unlabeled_category="Unknown"\n'
  r'\)\n'
  r'lvae\.train\(max_epochs=MAX_EPOCHS_SCANVI, batch_size=BATCH_SIZE, n_samples_per_label=2000\)\n'
  r'print\("  Training complete"\)\n\n'
  r'lvae\.save\(OUTPUT_DIR / "scanvi_model", overwrite=True\)\n'
  r'print\("  Model saved"\)',
  re.MULTILINE,
)
DOTPLOT_BLOCK_PATTERN = re.compile(
  r'dp = sc\.pl\.dotplot\(\n'
  r'\s+adata, var_names=available_flat, groupby=CELLTYPE_L3_REFINED_COL,\n'
  r"\s+standard_scale='var', use_raw=True, show=False, figsize=\(fig_w, fig_h\)\n"
  r'\)\n\n'
  r"dp\.fig\.suptitle\('Myeloid Validation - Refined Annotation', fontsize=18, y=0\.998, weight='bold'\)\n"
  r'dp\.fig\.tight_layout\(\)\n\n'
  r"dp\.savefig\(OUTPUT_DIR / 'dotplot_VALIDATION\.pdf', bbox_inches='tight'\)\n"
  r"dp\.savefig\(OUTPUT_DIR / 'dotplot_VALIDATION\.png', dpi=300, bbox_inches='tight'\)\n\n"
  r'print\("  Dotplot saved"\)',
  re.MULTILINE,
)
NEW_SCVI_BLOCK = '''print("Training scVI...")
scvi_model_path = OUTPUT_DIR / "scvi_model"
if scvi_model_path.exists():
  print("  Found existing scVI model; loading instead of retraining")
  vae = scvi.model.SCVI.load(scvi_model_path, adata=adata_scvi)
  print("  Loaded existing model")
else:
  vae.train(max_epochs=MAX_EPOCHS_SCVI, batch_size=BATCH_SIZE, early_stopping=True)
  print("  Training complete")
  vae.save(scvi_model_path, overwrite=True)
  print("  Model saved")'''
NEW_SCANVI_BLOCK = '''print("Training scANVI...")
scanvi_model_path = OUTPUT_DIR / "scanvi_model"
if scanvi_model_path.exists():
  print("  Found existing scANVI model; loading instead of retraining")
  lvae = scvi.model.SCANVI.load(scanvi_model_path, adata=adata_scvi)
  print("  Loaded existing model")
else:
  lvae = scvi.model.SCANVI.from_scvi_model(
    vae, adata=adata_scvi, labels_key=CELLTYPE_L3_REFINED_COL, unlabeled_category="Unknown"
  )
  lvae.train(
    max_epochs=MAX_EPOCHS_SCANVI,
    batch_size=BATCH_SIZE,
    n_samples_per_label=2000,
    early_stopping=True,
    early_stopping_patience=30,
  )
  print("  Training complete")
  lvae.save(scanvi_model_path, overwrite=True)
  print("  Model saved")'''
NEW_DOTPLOT_BLOCK = '''dp = sc.pl.dotplot(
    adata, var_names=available_flat, groupby=CELLTYPE_L3_REFINED_COL,
    standard_scale='var', use_raw=True, show=False, return_fig=True, figsize=(fig_w, fig_h)
)

if hasattr(dp, "make_figure"):
  dp.make_figure()

fig = getattr(dp, "fig", None)
if fig is None and isinstance(dp, dict):
  for value in dp.values():
    if hasattr(value, "figure"):
      fig = value.figure
      break

if fig is None:
  raise RuntimeError(f"Unexpected dotplot return type: {type(dp)!r}")

fig.suptitle('Myeloid Validation - Refined Annotation', fontsize=18, y=0.998, weight='bold')
fig.tight_layout()

fig.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.pdf', bbox_inches='tight')
fig.savefig(OUTPUT_DIR / 'dotplot_VALIDATION.png', dpi=300, bbox_inches='tight')
plt.close(fig)

print("  Dotplot saved")'''


def build_wrapped_source() -> str:
  if not BASE_SCRIPT.exists():
    raise FileNotFoundError(f"Base script not found: {BASE_SCRIPT}")

  source = BASE_SCRIPT.read_text(encoding="utf-8")
  if OLD_OUTPUT_LINE not in source:
    raise RuntimeError("Could not find the expected OUTPUT_DIR line in the base script")
  source = source.replace(OLD_OUTPUT_LINE, NEW_OUTPUT_LINE, 1)

  source, scvi_subs = SCVI_BLOCK_PATTERN.subn(NEW_SCVI_BLOCK, source, count=1)
  if scvi_subs != 1:
    raise RuntimeError("Could not find the expected scVI training block in the base script")

  source, scanvi_subs = SCANVI_BLOCK_PATTERN.subn(NEW_SCANVI_BLOCK, source, count=1)
  if scanvi_subs != 1:
    raise RuntimeError("Could not find the expected scANVI training block in the base script")

  source, dotplot_subs = DOTPLOT_BLOCK_PATTERN.subn(NEW_DOTPLOT_BLOCK, source, count=1)
  if dotplot_subs != 1:
    raise RuntimeError("Could not find the expected dotplot block in the base script")

  return source


def main() -> None:
  # Force CPU mode before the wrapped script imports torch.
  os.environ["CUDA_VISIBLE_DEVICES"] = ""
  os.environ.setdefault("JAX_PLATFORMS", "cpu")
  os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

  print("[wrapper] Running myeloid scVI/scANVI in CPU-safe mode (2026-04-16 rerun)")
  print("[wrapper] CUDA visibility disabled because the installed GPU driver is too old for the current PyTorch build")

  source = build_wrapped_source()
  exec_globals = {
    "__name__": "__main__",
    "__file__": str(BASE_SCRIPT),
  }
  exec(compile(source, str(BASE_SCRIPT), "exec"), exec_globals)


if __name__ == "__main__":
  main()
