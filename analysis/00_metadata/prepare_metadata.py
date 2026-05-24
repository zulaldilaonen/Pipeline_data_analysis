#!/usr/bin/env python3

"""
prepare_metadata.py

Purpose:
    Prepare a cleaned metadata table for downstream RNA-seq analysis.

How to run:
    cd ~/Pipeline_data_analysis/analysis
    python3 00_metadata/prepare_metadata.py
"""

import pandas as pd
from pathlib import Path

input_file = Path("../clean_csv_file/SraRunTable_clean.csv")
output_dir = Path("00_metadata")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file)

required_columns = ["Run", "AGE", "sex", "pmi", "brain_region", "diagnosis"]
missing_columns = [col for col in required_columns if col not in df.columns]
if missing_columns:
    raise ValueError(f"Missing required columns: {missing_columns}")

metadata = df[required_columns].copy()

metadata = metadata.rename(columns={
    "Run": "sample",
    "AGE": "age"
})

for col in metadata.columns:
    if metadata[col].dtype == object:
        metadata[col] = metadata[col].astype(str).str.strip()

metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce")
metadata["pmi"] = pd.to_numeric(metadata["pmi"], errors="coerce")

metadata["sex"] = metadata["sex"].astype(str).str.lower().str.strip()
metadata["brain_region"] = metadata["brain_region"].astype(str).str.strip()
metadata["diagnosis"] = metadata["diagnosis"].astype(str).str.strip()

metadata.to_csv(output_dir / "metadata_raw_selected.csv", index=False)

duplicated_samples = metadata["sample"][metadata["sample"].duplicated()].tolist()
metadata = metadata.drop_duplicates(subset="sample").copy()

with open(output_dir / "metadata_qc.txt", "w", encoding="utf-8") as f:
    f.write("Metadata QC Summary\n")
    f.write("===================\n\n")
    f.write(f"Total rows in original table: {len(df)}\n")
    f.write(f"Total unique samples: {len(metadata)}\n")
    f.write(f"Number of duplicated samples: {len(duplicated_samples)}\n\n")

    if duplicated_samples:
        f.write("Duplicated sample IDs:\n")
        for sample in duplicated_samples:
            f.write(f"- {sample}\n")
        f.write("\n")

    f.write("Missing values per column:\n")
    f.write(metadata.isna().sum().to_string())
    f.write("\n\n")

    for col in ["sex", "brain_region", "diagnosis"]:
        f.write(f"Value counts for {col}:\n")
        f.write(metadata[col].value_counts(dropna=False).to_string())
        f.write("\n\n")

metadata.to_csv(output_dir / "metadata.csv", index=False)

print("Files created successfully:")
print(output_dir / "metadata_raw_selected.csv")
print(output_dir / "metadata.csv")
print(output_dir / "metadata_qc.txt")
