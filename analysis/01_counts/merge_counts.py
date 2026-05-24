#!/usr/bin/env python3

"""
merge_counts.py

Purpose:
    Merge featureCounts output files into a single count matrix.

How to run:
    cd ~/Pipeline_data_analysis/analysis
    python 01_counts/merge_counts.py
"""

from pathlib import Path
import pandas as pd

# Resolve project paths safely
script_dir = Path(__file__).resolve().parent
analysis_dir = script_dir.parent
project_root = analysis_dir.parent
counts_dir = project_root / "counts"

output_matrix = script_dir / "merged_counts.tsv"
output_qc = script_dir / "count_merge_qc.txt"

# Find all count files
count_files = sorted(counts_dir.glob("*_counts.txt"))

if not count_files:
    raise FileNotFoundError(f"No *_counts.txt files found in {counts_dir}")

merged_df = None
sample_names = []

for file in count_files:
    sample_name = file.name.replace("_counts.txt", "")
    sample_names.append(sample_name)

    # Read featureCounts file (skip comment lines)
    df = pd.read_csv(file, sep="\t", comment="#")

    if "Geneid" not in df.columns:
        raise ValueError(f"'Geneid' column not found in {file}")

    # Last column is usually count column
    count_col = df.columns[-1]

    temp_df = df[["Geneid", count_col]].copy()
    temp_df.columns = ["Geneid", sample_name]

    if merged_df is None:
        merged_df = temp_df
    else:
        merged_df = merged_df.merge(temp_df, on="Geneid", how="outer")

# Fill missing values with 0
merged_df = merged_df.fillna(0)

# Convert all counts to integer
for col in merged_df.columns[1:]:
    merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce").fillna(0).astype(int)

# Save merged matrix
merged_df.to_csv(output_matrix, sep="\t", index=False)

# Write QC summary
with open(output_qc, "w") as f:
    f.write("Count Merge QC Summary\n")
    f.write("======================\n\n")
    f.write(f"Counts directory: {counts_dir}\n")
    f.write(f"Number of count files: {len(count_files)}\n")
    f.write(f"Number of genes: {merged_df.shape[0]}\n")
    f.write(f"Number of samples: {merged_df.shape[1] - 1}\n\n")

    f.write("Sample names:\n")
    for s in sample_names:
        f.write(f"- {s}\n")

print("Merge completed successfully")
print(f"Output: {output_matrix}")
print(f"QC: {output_qc}")
