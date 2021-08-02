import pandas as pd
import numpy as np


def clean_gene_df(df):
    ## if cadd is missing (nan) give it the mean of the same type (Delet, Mis)
    mean_cadd_delet = np.mean(
        df.loc[(df["Deleterious"] == 1) & (~df["CADD"].isna())]["CADD"].values
    )

    mean_cadd_missense = np.mean(
        df.loc[(df["Missense"] == 1) & (~df["CADD"].isna())]["CADD"].values
    )

    mean_cadd_synon = np.mean(
        df.loc[(df["Synonymous"] == 1) & (~df["CADD"].isna())]["CADD"].values
    )

    if np.isnan(mean_cadd_synon) or pd.isnull(mean_cadd_synon):
        mean_cadd_synon = 0

    mean_gerp = np.mean(
        df.loc[
            (
                (df["Deleterious"] == 1) |
                (df["Missense"] == 1) |
                (df["Synonymous"] == 1)
            ) &
            (~df["GERP"].isna())
        ]["GERP"].values
    )

    mean_phylop = np.mean(
        df.loc[
            ((df["Deleterious"] == 1) | (df["Missense"] == 1)) &
            (~df["phyloP"].isna())
        ]["phyloP"].values
    )

    cleaned_cadd = []
    cleaned_gerp = []
    cleaned_phylop = []
    cleaned_af = []
    for index, row in df.iterrows():
        if pd.isnull(row["CADD"]):
            if row["Missense"] == 1:
                cleaned_cadd.append(mean_cadd_missense)
            elif row["Deleterious"] == 1:
                cleaned_cadd.append(mean_cadd_delet)
            elif row["Synonymous"] == 1:
                cleaned_cadd.append(mean_cadd_synon)
        else:
            cleaned_cadd.append(row["CADD"])

        if pd.isnull(row["GERP"]):
            cleaned_gerp.append(mean_gerp)
        else:
            cleaned_gerp.append(row["GERP"])

        if pd.isnull(row["phyloP"]):
            cleaned_phylop.append(mean_phylop)
        else:
            cleaned_phylop.append(row["phyloP"])

        if pd.isnull(row["allele_frequency"]):
            cleaned_af.append(-6)
        else:
            cleaned_af.append(row["allele_frequency"])


    df["CADD"] = cleaned_cadd
    df["GERP"] = cleaned_gerp
    df["phyloP"] = cleaned_phylop
    df["allele_frequency"] = cleaned_af

    return df
    # df["CADD"] = df["CADD"].apply(lambda x : )

if __name__ == "__main__":
    genes = ["BRCA2", "BRCA1", "LDLR", "APOB", "PCSK9"]
    for g in genes:
        df = pd.read_csv(f"{g}.csv")
        clean_gene_df(df)
        df.to_csv(f"cleaned_{g}.csv")

        # df.to_csv("cleaned_BRCA2.csv")
