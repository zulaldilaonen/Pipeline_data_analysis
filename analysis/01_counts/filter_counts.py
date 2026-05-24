#!/usr/bin/env python3

"""
filter_counts.py

Purpose:
    Filter merged count matrix based on kept samples.

How to run:
    python 01_counts/filter_counts.py
"""

import pandas as pd

counts = pd.read_csv("01_counts/merged_counts.tsv", sep="\t")
kept = pd.read_csv("02_qc/kept_samples_v3.tsv", sep="\t")

cols_to_keep = ["Geneid"] + kept["sample"].tolist()

filtered = counts[cols_to_keep]

filtered.to_csv("01_counts/merged_counts_filtered.tsv", sep="\t", index=False)

print("merged_counts_filtered.tsv created")
print("Genes:", filtered.shape[0])
print("Samples:", filtered.shape[1] - 1)
