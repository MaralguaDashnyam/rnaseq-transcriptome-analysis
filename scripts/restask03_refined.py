# -*- coding: utf-8 -*-
"""
ResTask03 — Gene/Transcript Comparison (Reference vs StringTie) — Refined
Author: Hera Dashnyam (with assist)
Purpose:
  - Produce CLEAN, comparable tables for reference genes and StringTie results.
  - Add strand + category cols; normalize IDs (remove "gene:" prefixes).
  - Add transcript-level categories; rank transcripts by length (short → long).
  - Generate gene-level summary for StringTie (C1/C2/C3) using overlaps.
  - Export a unified comparison table with same columns for quick diffing.

How to run (example — edit paths to your files):
    python restask03_refined.py \
        --assembled /path/to/YapPool_merged.gtf \
        --reference /path/to/Gallus_gallus.GRCg7b.113.gff3 \
        --outdir /path/to/ResTask03_out

Outputs (CSV):
  1) reference_genes.csv
  2) assembled_transcripts.csv
  3) stringtie_gene_summary.csv
  4) gene_comparison_unified.csv   # SAME columns for ref vs stringtie
  5) loci_C1.bed / loci_C2.bed / loci_C3.bed

Notes:
  - Overlaps are computed strand-agnostic by default (can toggle --strand_aware).
  - "category" in transcript table: NOVEL (no overlap) or OVERLAP.
  - "category" in gene tables: C1 (purely novel locus), C2, C3 (see below).
  - Dominant strand for a StringTie locus = majority transcript strand; ties=".".
"""

import argparse, gzip, sys, os
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

# Optional: pybedtools (bedtools must be available in PATH)
try:
    import pybedtools
except Exception as e:
    print("ERROR: pybedtools is required. pip install pybedtools", file=sys.stderr)
    raise

pd.set_option("display.max_columns", 200)
pd.set_option("display.width", 180)


# ---------------------------
# Helpers
# ---------------------------
def _open_auto(p: Path):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

def _parse_attrs(attr: str):
    d = {}
    for part in str(attr).strip().split(";"):
        part = part.strip()
        if not part: 
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            continue
        d[k] = v.strip().strip('"')
    return d

def _coerce_int(x):
    try: 
        return int(x)
    except Exception:
        return np.nan

def sanitize_intervals(df, start_col, end_col, chrom_col="chrom", name_cols=None, label=""):
    d = df.copy()
    before = len(d)

    d[chrom_col] = d[chrom_col].astype(str).str.strip()
    d[start_col] = d[start_col].map(_coerce_int)
    d[end_col]   = d[end_col].map(_coerce_int)

    d = d.dropna(subset=[chrom_col, start_col, end_col])
    d[start_col] = d[start_col].astype(int)
    d[end_col]   = d[end_col].astype(int)

    bad_end = d[end_col] < d[start_col]
    d = d[~bad_end]

    if name_cols:
        for c in name_cols:
            if c in d.columns:
                d[c] = d[c].fillna(".").astype(str)

    print(f"[sanitize {label}] in={before} kept={len(d)} dropped={before-len(d)}")
    return d.reset_index(drop=True)

def to_bed(df, path, name_cols=("gene_id","n_transcripts")):
    with open(path, "w") as fh:
        for _, r in df.iterrows():
            chrom = str(r["chrom"])
            start0 = max(int(r["start"]) - 1, 0)
            end1   = int(r["end"])
            label  = "|".join(f"{k}={r[k]}" for k in name_cols if k in r and pd.notna(r[k]))
            if not label:
                label = str(r.get("gene_id","."))
            fh.write(f"{chrom}\t{start0}\t{end1}\t{label}\n")


# ---------------------------
# Parsing: Reference (GFF3/GTF)
# ---------------------------
def parse_reference_genes(reference_path: Path):
    """Return gene table with: ref_gene_id, chrom, strand, ref_start, ref_end, ref_span_bp, ref_n_transcripts, category='REFERENCE'"""
    genes = {}  # gid -> (chrom, strand, min_start, max_end)
    tx_count = defaultdict(int)

    with _open_auto(reference_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, src, ftype, start, end, score, strand, phase, attr = cols
            start, end = int(start), int(end)
            a = _parse_attrs(attr)

            if ftype.lower() in ("gene","pseudogene"):
                gid = a.get("ID") or a.get("gene_id") or a.get("Name")
                if not gid:
                    continue
                if gid in genes:
                    c, s, s0, e0 = genes[gid]
                    genes[gid] = (chrom, strand if strand != "." else s, min(s0, start), max(e0, end))
                else:
                    genes[gid] = (chrom, strand, start, end)

            elif ftype.lower() in ("mrna","transcript"):
                gid = a.get("Parent") or a.get("gene_id") or a.get("gene") or a.get("gene_name")
                if gid:
                    tx_count[gid] += 1

    rows = []
    for gid, (chrom, strand, s, e) in genes.items():
        clean_gid = gid.replace("gene:", "")
        rows.append([clean_gid, chrom, strand, s, e, int(e - s + 1), int(tx_count.get(gid, 0)), "REFERENCE"])

    ref_df = pd.DataFrame(rows, columns=[
        "gene_id","chrom","strand","start","end","span_bp","n_transcripts","category"
    ]).sort_values(["chrom","start","end"]).reset_index(drop=True)

    ref_df = sanitize_intervals(ref_df, "start","end","chrom", name_cols=["gene_id","strand","category"], label="ref_df")
    return ref_df


# ---------------------------
# Parsing: Assembled GTF (StringTie)
# ---------------------------
def parse_assembled_gtf(gtf_path: Path):
    """Return transcript table: chrom, start, end, tx_id, strand, assembled_gene_id"""
    bt = pybedtools.BedTool(str(gtf_path))

    def parse_attrs_gtf(field):
        d = {}
        for part in str(field).strip().split(";"):
            part = part.strip()
            if not part:
                continue
            if " " in part and "=" not in part:
                k, v = part.split(" ", 1); d[k] = v.strip('"; ')
            elif "=" in part:
                k, v = part.split("=", 1); d[k] = v.strip('"; ')
        return d

    tx_rows, exon_rows = [], []
    for f in bt:
        ftype = f[2]
        attrs = parse_attrs_gtf(f[8] if len(f) > 8 else "")
        tx_id = attrs.get("transcript_id") or attrs.get("ID") or attrs.get("transcript") or attrs.get("Name") or attrs.get("Parent")
        gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("gene_name") or attrs.get("geneID") or attrs.get("Parent")
        if ftype in ("transcript","mRNA"):
            tx_rows.append([f.chrom, int(f.start), int(f.end), tx_id or ".", f.strand, gene_id or "."])
        elif ftype == "exon":
            exon_rows.append([f.chrom, int(f.start), int(f.end), tx_id or ".", f.strand, gene_id or "."])

    tx_df = pd.DataFrame(tx_rows, columns=["chrom","start","end","tx_id","strand","assembled_gene_id"])
    exon_df = pd.DataFrame(exon_rows, columns=["chrom","start","end","tx_id","strand","assembled_gene_id"])

    # Fallback: infer transcript bounds from exons if transcript rows are missing
    if tx_df.empty and not exon_df.empty:
        agg = exon_df.groupby(["tx_id","assembled_gene_id","chrom","strand"], dropna=False)\
                     .agg(start=("start","min"), end=("end","max")).reset_index()
        tx_df = agg[["chrom","start","end","tx_id","strand","assembled_gene_id"]].copy()

    tx_df = sanitize_intervals(
        tx_df, start_col="start", end_col="end", chrom_col="chrom",
        name_cols=["tx_id","assembled_gene_id","strand"], label="tx_df"
    )
    tx_df["tx_len"] = (tx_df["end"] - tx_df["start"] + 1).astype(int)
    return tx_df, exon_df


# ---------------------------
# Overlap logic
# ---------------------------
def build_bed(df, name_col, start_col="start", end_col="end"):
    ivals = []
    for _, r in df.sort_values(["chrom", start_col, end_col]).iterrows():
        start0 = max(int(r[start_col]) - 1, 0)
        end1   = int(r[end_col])
        name   = str(r[name_col]) if pd.notna(r[name_col]) else "."
        ivals.append(pybedtools.Interval(str(r["chrom"]), start0, end1, name=name))
    return pybedtools.BedTool(ivals).sort()

def compute_overlaps(tx_df, ref_df, strand_aware=False):
    asm_bt = build_bed(tx_df, "tx_id")
    ref_bt = build_bed(ref_df.rename(columns={"gene_id":"ref_gene_id"}), "ref_gene_id")

    inter = asm_bt.intersect(ref_bt, wa=True, wb=True, s=strand_aware)
    tx2ref = {}
    for iv in inter:
        tx_name = iv[3]
        ref_gid = iv[7]
        tx2ref.setdefault(tx_name, set()).add(ref_gid)

    rows = []
    for _, r in tx_df.sort_values(["chrom","start","end"]).iterrows():
        txid = r["tx_id"]
        ov = tx2ref.get(txid, set())
        cat = "NOVEL" if len(ov) == 0 else "OVERLAP"
        rows.append({
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "tx_id": txid,
            "strand": r["strand"],
            "assembled_gene_id": r["assembled_gene_id"],
            "tx_len": int(r["tx_len"]),
            "overlap_ref_gene_ids": ";".join(sorted(ov)) if ov else "",
            "category": cat,
        })
    out = pd.DataFrame(rows)
    out = out.sort_values(["chrom","assembled_gene_id","tx_len","start","end"]).reset_index(drop=True)
    return out


# ---------------------------
# Cluster transcripts into loci (StringTie "gene-level")
# ---------------------------
def cluster_stringtie_loci(tx_table, strand_aware=False):
    asm_bt = build_bed(tx_table.rename(columns={"tx_id":"name"}).assign(name=lambda d: d["tx_id"]), "name")
    clustered = asm_bt.cluster(s=strand_aware)

    # tx_id -> cluster_id
    cluster_map = {iv[3]: iv[-1] for iv in clustered}
    tx_table = tx_table.copy()
    tx_table["cluster_id"] = tx_table["tx_id"].map(cluster_map)

    # gene-level aggregation (per cluster)
    loci = []
    for cid, grp in tx_table.groupby("cluster_id"):
        chrom = grp["chrom"].iloc[0]
        start = int(grp["start"].min())
        end   = int(grp["end"].max())
        n_tx  = int(grp.shape[0])
        # dominant strand
        c = Counter(grp["strand"].astype(str))
        if not c:
            dom = "."
        else:
            top = c.most_common(2)
            dom = top[0][0] if len(top)==1 or top[0][1] > top[1][1] else "."

        # any reference overlaps?
        overlaps = set()
        for s in grp.get("overlap_ref_gene_ids",""):
            if isinstance(s,str) and s:
                overlaps.update(s.split(";"))
        overlaps = {g for g in overlaps if g}

        # Final gene-level category:
        # C1: all transcripts novel (no overlap)
        # C2: overlaps reference, and assembled gene_id does NOT match the overlapped ref gene IDs
        # C3: overlaps reference, and assembled gene_id matches one of the ref gene IDs
        cats = set(grp["category"])
        if cats == {"NOVEL"} and len(overlaps) == 0:
            final_cat = "C1"
        else:
            asm_gene_ids = set(grp["assembled_gene_id"].astype(str))
            if len(overlaps) >= 1:
                final_cat = "C3" if len(asm_gene_ids.intersection(overlaps)) >= 1 else "C2"
            else:
                final_cat = "C1"

        loci.append({
            "gene_id": f"ASMLOC_{cid}",
            "chrom": chrom,
            "strand": dom,
            "start": start,
            "end": end,
            "span_bp": int(end - start + 1),
            "n_transcripts": n_tx,
            "overlap_ref_gene_ids": ";".join(sorted(overlaps)),
            "category": final_cat
        })

    gene_df = pd.DataFrame(loci).sort_values(["chrom","start","end"]).reset_index(drop=True)
    return gene_df


# ---------------------------
# Unified comparison table
# ---------------------------
def build_unified(ref_genes, asm_genes):
    ref_norm = ref_genes.rename(columns={"gene_id":"gene_id","span_bp":"span_bp","n_transcripts":"n_transcripts"}).copy()
    ref_norm["source"] = "REFERENCE"

    asm_norm = asm_genes.rename(columns={"gene_id":"gene_id","span_bp":"span_bp","n_transcripts":"n_transcripts"}).copy()
    asm_norm["source"] = "STRINGTIE"

    cols = ["source","gene_id","chrom","strand","start","end","span_bp","n_transcripts","category"]
    uni = pd.concat([ref_norm[cols], asm_norm[cols]], ignore_index=True)\
            .sort_values(["chrom","start","end","source"]).reset_index(drop=True)
    return uni


# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembled", required=True, type=Path, help="StringTie merged GTF")
    ap.add_argument("--reference", required=True, type=Path, help="Reference GFF3/GTF")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory")
    ap.add_argument("--strand_aware", action="store_true", help="Use strand-aware overlap/cluster")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # 1) Reference genes
    ref_genes = parse_reference_genes(args.reference)

    # 2) Assembled transcripts
    tx_df, _ = parse_assembled_gtf(args.assembled)

    # 3) Transcript vs reference overlaps + categories
    tx_cat = compute_overlaps(tx_df, ref_genes, strand_aware=args.strand_aware)

    # 4) StringTie gene-level loci & categories
    asm_genes = cluster_stringtie_loci(tx_cat, strand_aware=args.strand_aware)

    # 5) Unified gene table for side-by-side comparison
    unified = build_unified(ref_genes, asm_genes)

    # 6) Save CSVs
    ref_csv = args.outdir / "reference_genes.csv"
    tx_csv  = args.outdir / "assembled_transcripts.csv"
    asm_csv = args.outdir / "stringtie_gene_summary.csv"
    uni_csv = args.outdir / "gene_comparison_unified.csv"

    ref_genes.to_csv(ref_csv, index=False)
    tx_cat.sort_values(["assembled_gene_id","tx_len","start"]).to_csv(tx_csv, index=False)
    asm_genes.to_csv(asm_csv, index=False)
    unified.to_csv(uni_csv, index=False)

    print(f"Saved: {ref_csv}")
    print(f"Saved: {tx_csv}")
    print(f"Saved: {asm_csv}")
    print(f"Saved: {uni_csv}")

    # 7) IGV BED files per StringTie locus category
    to_bed(asm_genes[asm_genes["category"]=="C1"][["chrom","start","end","gene_id","n_transcripts"]].rename(columns={"gene_id":"gene_id"}),
           args.outdir / "loci_C1.bed",
           name_cols=("gene_id","n_transcripts"))
    to_bed(asm_genes[asm_genes["category"]=="C2"][["chrom","start","end","gene_id","n_transcripts"]].rename(columns={"gene_id":"gene_id"}),
           args.outdir / "loci_C2.bed",
           name_cols=("gene_id","n_transcripts"))
    to_bed(asm_genes[asm_genes["category"]=="C3"][["chrom","start","end","gene_id","n_transcripts"]].rename(columns={"gene_id":"gene_id"}),
           args.outdir / "loci_C3.bed",
           name_cols=("gene_id","n_transcripts"))

    print("Saved BEDs: loci_C1.bed / loci_C2.bed / loci_C3.bed")

if __name__ == "__main__":
    main()
