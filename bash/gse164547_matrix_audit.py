#!/usr/bin/env python3
import os, re, sys, gzip, itertools, collections, argparse
from datetime import datetime

def split_fields(s: str, mode: str):
    if mode == "tab":
        a = s.split("\t")
    elif mode == "comma":
        a = s.split(",")
    elif mode == "ws":
        a = re.split(r"\s+", s.strip())
    else:
        raise ValueError(mode)
    # drop trailing empty fields (common if line ends with delimiter)
    while a and a[-1] == "":
        a.pop()
    return a

def detect_delim(b: bytes):
    # Pick the delimiter with the highest count in line1 among TAB, comma, space.
    # If TAB exists at all, usually safest for these matrices.
    tab = b.count(b"\t")
    com = b.count(b",")
    sp  = b.count(b" ")
    nul = b.count(b"\x00")

    if tab > 0 and tab >= com and tab >= sp:
        mode = "tab"
    elif com > 0 and com >= sp:
        mode = "comma"
    else:
        mode = "ws"
    return mode, tab, com, sp, nul

def read_first_n_lines_gz(path, n=3):
    with gzip.open(path, "rb") as fh:
        lines = [fh.readline() for _ in range(n)]
    return lines

def visible_preview(b: bytes, limit=240):
    s = b.decode("utf-8", "replace")
    s = s.replace("\t","<TAB>").replace("\r","<CR>").replace("\n","<LF>")
    return s[:limit] + ("..." if len(s) > limit else "")

def audit_one_file(path, n_lines=200):
    # read raw first 3 lines
    b1, b2, b3 = read_first_n_lines_gz(path, 3)

    mode, tab, com, sp, nul = detect_delim(b1)
    l1 = b1.decode("utf-8","replace").rstrip("\n")
    l2 = b2.decode("utf-8","replace").rstrip("\n")

    f1 = split_fields(l1, mode)
    f2 = split_fields(l2, mode)

    # field count distribution in first n_lines
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        lines = list(itertools.islice(fh, n_lines))
    dist = collections.Counter(len(split_fields(line.rstrip("\n"), mode)) for line in lines)
    top_dist = dist.most_common(5)

    # heuristics
    issues = []
    reason = "ok"

    # NUL byte indicates binary corruption or wrong decompression at some point
    if nul > 0:
        issues.append("contains_NUL_bytes")
        reason = "contains_NUL_bytes"

    # column consistency: if >1 distinct field counts, it's inconsistent
    if len(dist) > 1:
        issues.append("inconsistent_field_count_first200")
        if reason == "ok":
            reason = "matrix_cols_inconsistent_head"

    # header missing gene placeholder: common pattern L2 has exactly +1 field than L1
    # and L2[0] looks like a gene symbol (non-numeric, often contains letters) while L1[0] looks like a barcode
    # We only flag if dist is otherwise consistent (or at least dominant counts) to avoid false positives.
    if len(f1) > 1 and len(f2) > 1 and (len(f2) - len(f1) == 1):
        issues.append("header_missing_gene_placeholder")
        if reason == "ok":
            reason = "header_missing_gene_placeholder"

    # trailing delimiter: if raw line ends with delimiter -> produces empty last field before trimming
    raw_line1 = b1.decode("utf-8","replace").rstrip("\n")
    if mode == "tab" and raw_line1.endswith("\t"):
        issues.append("line1_trailing_tab")
        if reason == "ok":
            reason = "line_trailing_delim"
    if mode == "comma" and raw_line1.endswith(","):
        issues.append("line1_trailing_comma")
        if reason == "ok":
            reason = "line_trailing_delim"

    return {
        "path": path,
        "mode": mode,
        "tab": tab,
        "comma": com,
        "space": sp,
        "nul": nul,
        "n_field_L1": len(f1),
        "n_field_L2": len(f2),
        "diff_L2_L1": len(f2) - len(f1),
        "dist_top5": top_dist,
        "dist_unique": len(dist),
        "reason": reason,
        "issues": issues,
        "preview_L1": visible_preview(b1),
        "preview_L2": visible_preview(b2),
        "preview_L3": visible_preview(b3),
    }

def iter_gz_files(root, dataset_id):
    base = os.path.join(root, dataset_id)
    if not os.path.isdir(base):
        return
    for gsm in sorted(os.listdir(base)):
        if not re.match(r"^(GSM|HRR)\d+$", gsm):
            continue
        d = os.path.join(base, gsm)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".gz") and re.search(r"(expression_matrix|counts|matrix)\.(txt|tsv|csv)\.gz$", fn, flags=re.I):
                yield dataset_id, gsm, os.path.join(d, fn)

def write_report(results, out_txt, out_tsv):
    # TXT report
    with open(out_txt, "w", encoding="utf-8") as out:
        out.write("============================================================\n")
        out.write("GSE164547 Matrix Audit Report\n")
        out.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")
        out.write("============================================================\n\n")

        for r in results:
            out.write("------------------------------------------------------------\n")
            out.write(f"DATASET: {r['dataset_id']}  SAMPLE: {r['sample_id']}\n")
            out.write(f"FILE: {r['path']}\n")
            out.write(f"DELIM_MODE: {r['mode']}  tabs={r['tab']} commas={r['comma']} spaces={r['space']} nul={r['nul']}\n")
            out.write(f"FIELDS: L1={r['n_field_L1']}  L2={r['n_field_L2']}  diff(L2-L1)={r['diff_L2_L1']}\n")
            out.write(f"DIST(first200): unique={r['dist_unique']} top5={r['dist_top5']}\n")
            out.write(f"REASON: {r['reason']}\n")
            if r["issues"]:
                out.write("ISSUES: " + ", ".join(r["issues"]) + "\n")
            out.write("\n")
            out.write("PREVIEW L1: " + r["preview_L1"] + "\n")
            out.write("PREVIEW L2: " + r["preview_L2"] + "\n")
            out.write("PREVIEW L3: " + r["preview_L3"] + "\n\n")

    # TSV summary
    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("\t".join([
            "dataset_id","sample_id","file","delim_mode","nul_bytes",
            "fields_L1","fields_L2","diff_L2_L1",
            "dist_unique_first200","dist_top5_first200",
            "reason","issues"
        ]) + "\n")
        for r in results:
            out.write("\t".join([
                r["dataset_id"],
                r["sample_id"],
                r["path"],
                r["mode"],
                str(r["nul"]),
                str(r["n_field_L1"]),
                str(r["n_field_L2"]),
                str(r["diff_L2_L1"]),
                str(r["dist_unique"]),
                str(r["dist_top5"]),
                r["reason"],
                ",".join(r["issues"]) if r["issues"] else ""
            ]) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_root", required=True)
    ap.add_argument("--dataset_id", default="GSE164547")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n_lines", type=int, default=200)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    out_txt = os.path.join(args.outdir, f"{args.dataset_id}.matrix_audit.{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
    out_tsv = os.path.join(args.outdir, f"{args.dataset_id}.matrix_audit_summary.{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv")

    results = []
    for dataset_id, sample_id, fpath in iter_gz_files(args.data_root, args.dataset_id):
        r = audit_one_file(fpath, n_lines=args.n_lines)
        r["dataset_id"] = dataset_id
        r["sample_id"] = sample_id
        results.append(r)

    if not results:
        print("No matching *.txt.gz/*.tsv.gz/*.csv.gz matrices found under:", os.path.join(args.data_root, args.dataset_id))
        sys.exit(2)

    # Sort: problematic first
    rank = {"ok": 9}
    def keyfun(r):
        return (0 if r["reason"] != "ok" else 1, r["reason"], r["sample_id"])
    results.sort(key=keyfun)

    write_report(results, out_txt, out_tsv)

    # Console summary
    from collections import Counter
    c = Counter(r["reason"] for r in results)
    print("============================================================")
    print("DONE")
    print("OUT TXT:", out_txt)
    print("OUT TSV:", out_tsv)
    print("Reason distribution:")
    for k,v in c.most_common():
        print(f"  {k}\t{v}")
    print("============================================================")

if __name__ == "__main__":
    main()



