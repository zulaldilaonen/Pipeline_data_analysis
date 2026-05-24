#!/usr/bin/env python3

"""
filter_metadata.py

Purpose:
    Filter metadata based on kept sample list.

How to run:
    python 00_metadata/filter_metadata.py
"""

import pandas as pd

meta = pd.read_csv("00_metadata/metadata.csv")
kept = pd.read_csv("02_qc/kept_samples_v3.tsv", sep="\t")

filtered = meta[meta["sample"].isin(kept["sample"])]

filtered.to_csv("00_metadata/metadata_filtered.csv", index=False)

print("metadata_filtered.csv created")
print("Samples:", len(filtered))
