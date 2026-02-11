import pandas as pd
import re

df = pd.read_csv("fh_LDLR_clinvar_full.csv")


# REMOVE THE ONE NOT CHECKED

# Clinical label
def map_label(x):
    if pd.isna(x):
        return None
    x = x.lower()
    if "pathogenic" in x:
        return 1
    if "benign" in x:
        return 0
    return None


df["y"] = df["ClinicalSignificance"].apply(map_label)
df = df[df["y"].notna()].copy()


# Variant parsing
def extract_protein(name):
    m = re.search(r"\(p\.([^)]+)\)", str(name))
    return m.group(1) if m else None


def extract_cdna(name):
    m = re.search(r":c\.([^ ]+)", str(name))
    return m.group(1) if m else None


# remove somatic
df = df[df["OriginSimple"].isin(["germline", "germline/somatic", "unknown", "not applicable"])].copy()

df["Protein_variant"] = df["Name"].apply(extract_protein)
df["CDNA_variant"] = df["Name"].apply(extract_cdna)

# FH label
pl = df["PhenotypeList"].str.lower().fillna("")

df["FH_label"] = pl.str.contains("amilial", regex=True).astype(int)
df["FH_type"] = "non-FH"
df.loc[df["FH_label"] == 1, "FH_type"] = "FH_unspecified"
df.loc[pl.str.contains("homo"), "FH_type"] = "FH_homozygous"
df.loc[pl.str.contains("hetero"), "FH_type"] = "FH_heterozygous"

final_cols = [
    "AlleleID",
    "Protein_variant",
    "CDNA_variant",
    "Start",
    "Stop",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "Type",
    "FH_label",
    "FH_type",
    "y",
]

bio_df = df[final_cols].dropna(subset=["y"])

bio_df.to_csv("LDLR_big.csv", index=False)
