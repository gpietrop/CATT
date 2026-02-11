import pandas as pd
import re


type_encoding = {
    "single nucleotide variant": 0,
    "Deletion": 1,
    "Duplication": 2,
    "Indel": 3,
    "Insertion": 4,
    "Microsatellite": 5,
}


def generate_small_ds(df, FH_flag=0):
    # extract only digits from CDNA_variant
    df["CDNA_variant"] = (df["CDNA_variant"].astype(str).str.extract(r"(\d+)").astype(float))

    # encode type
    df["Type"] = df["Type"].map(type_encoding)

    # drop nan
    df = df.dropna(subset=["Protein_variant", "CDNA_variant"])

    if FH_flag:
        df = df[df["FH_label"] == 1].copy()

    return df


df = pd.read_csv("LDLR_big.csv")

small_df = generate_small_ds(df=df)
small_df.to_csv("LDLR_small.csv", index=False)
