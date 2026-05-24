#!/usr/bin/env python3

"""
summarize_featurecounts.py

Purpose:
    Summarize featureCounts assignment statistics and filter samples
    based on assigned percentage threshold.

How to run:
    cd ~/Pipeline_data_analysis/analysis
    python 02_qc/summarize_featurecounts.py
"""

from pathlib import Path
import pandas as pd

# Resolve paths
script_dir = Path(__file__).resolve().parent
analysis_dir = script_dir.parent
project_root = analysis_dir.parent
counts_dir = project_root / "counts"

output_qc = script_dir / "featurecounts_qc.tsv"
output_kept = script_dir / "kept_samples.tsv"
output_removed = script_dir / "removed_samples.tsv"
output_summary = script_dir / "qc_summary.txt"

threshold = 35.0

summary_files = sorted(counts_dir.glob("*_counts.txt.summary"))

if not summary_files:
    raise FileNotFoundError(f"No *_counts.txt.summary files found in {counts_dir}")

records = []

for file in summary_files:
    sample = file.name.replace("_counts.txt.summary", "")
    df = pd.read_csv(file, sep="\t", index_col=0)

    if df.shape[1] != 1:
        raise ValueError(f"Unexpected summary format in {file}")

    sample_col = df.columns[0]

    assigned_reads = df.loc["Assigned", sample_col] if "Assigned" in df.index else 0
    total_reads = df[sample_col].sum()
    assigned_percent = (assigned_reads / total_reads) * 100 if total_reads > 0 else 0

    records.append({
        "sample": sample,
        "assigned_reads": int(assigned_reads),
        "total_reads": int(total_reads),
        "assigned_percent": round(assigned_percent, 4),
        "status": "keep" if assigned_percent >= threshold else "remove"
    })

qc_df = pd.DataFrame(records).sort_values("assigned_percent")

qc_df.to_csv(output_qc, sep="\t", index=False)

kept_df = qc_df[qc_df["status"] == "keep"].copy()
removed_df = qc_df[qc_df["status"] == "remove"].copy()

kept_df[["sample"]].to_csv(output_kept, sep="\t", index=False)
removed_df[["sample"]].to_csv(output_removed, sep="\t", index=False)

with open(output_summary, "w", encoding="utf-8") as f:
    f.write("FeatureCounts QC Summary\n")
    f.write("========================\n\n")
    f.write(f"Counts directory: {counts_dir}\n")
    f.write(f"Number of summary files: {len(summary_files)}\n")
    f.write(f"Threshold (assigned percent): {threshold}%\n")
    f.write(f"Kept samples: {len(kept_df)}\n")
    f.write(f"Removed samples: {len(removed_df)}\n\n")

    if not removed_df.empty:
        f.write("Removed sample list:\n")
        for sample in removed_df["sample"]:
            f.write(f"- {sample}\n")

print("Files created successfully:")
print(output_qc)
print(output_kept)
print(output_removed)
print(output_summary)
