import pandas as pd

# Silinecek GSM'ler
gsm_to_remove = {
    "GSM5292838","GSM5292839","GSM5292840","GSM5292841","GSM5292842","GSM5292843",
    "GSM5292844","GSM5292845","GSM5292846","GSM5292847","GSM5292848","GSM5292849",
    "GSM5292850","GSM5292851","GSM5292852","GSM5292853","GSM5292854","GSM5292855",
    "GSM5292856","GSM5292857","GSM5292858","GSM5292859","GSM5292860","GSM5292861",
    "GSM5292862","GSM5292863","GSM5292864","GSM5292865","GSM5292866","GSM5292867",
    "GSM5292868","GSM5292869","GSM5292870","GSM5292871","GSM5292872","GSM5292873",
    "GSM5292874","GSM5292875","GSM5292876"
}

# CSV oku
df = pd.read_csv("SraRunTable.csv", dtype=str).fillna("")

# GSM kolonunu bul
gsm_col = None
for col in df.columns:
    if df[col].str.contains(r"^GSM\d+$", regex=True).any():
        gsm_col = col
        break

if gsm_col is None:
    raise ValueError(f"GSM kolonu bulunamadı. Kolonlar: {list(df.columns)}")

# Run/SRR kolonunu bul
run_col = None
candidate_run_cols = ["Run", "Run_s", "SRR", "Accession"]
for col in candidate_run_cols:
    if col in df.columns and df[col].str.contains(r"^SRR\d+$", regex=True).any():
        run_col = col
        break

if run_col is None:
    for col in df.columns:
        if df[col].str.contains(r"^SRR\d+$", regex=True).any():
            run_col = col
            break

if run_col is None:
    raise ValueError(f"SRR/Run kolonu bulunamadı. Kolonlar: {list(df.columns)}")

# Filtrele
removed_mask = df[gsm_col].isin(gsm_to_remove)
df_clean = df.loc[~removed_mask].copy()

# Temiz CSV yaz
df_clean.to_csv("SraRunTable_clean.csv", index=False)

# Temiz SRR listesi yaz
srr_list = (
    df_clean[run_col]
    .dropna()
    .astype(str)
    .str.strip()
)
srr_list = srr_list[srr_list.str.match(r"^SRR\d+$", na=False)]

with open("SRR_Acc_List_clean.txt", "w") as f:
    for srr in srr_list:
        f.write(srr + "\n")

print(f"GSM kolonu: {gsm_col}")
print(f"Run kolonu: {run_col}")
print(f"Silinen satır sayısı: {removed_mask.sum()}")
print(f"Kalan satır sayısı: {len(df_clean)}")
print("Oluşturulan dosyalar:")
print("- SraRunTable_clean.csv")
print("- SRR_Acc_List_clean.txt")