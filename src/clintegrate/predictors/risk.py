from clintegrate.predictors.constants import (
    included_genes,
    help_statement,
    gene_to_required_fields,
    default_variant_info_dictionaries
)
from clintegrate.predictors.Exceptions  import (
    InvalidFormatError,
    VariantFormatException,
    InvalidGeneException,
    NotInitializedException
)
from clintegrate.predictors.tools import (
    validate_variant_format
 )
import pandas as pd
import pickle
import pkg_resources


class IntegrativePredictiveModel():
    def __init__(self):
        self.gene = None
        self.model = None
        self.precomputed_snps = None
        self.initialized = False

    def load_df(self, p, use_index_0=False):
        stream = pkg_resources.resource_stream(__name__, p)
        if use_index_0:
            return pd.read_csv(stream, encoding='latin-1', index_col = 0)
        else:
            return pd.read_csv(stream, encoding='latin-1')

    def initialize(self, gene):
        if gene in included_genes:
            self.gene = gene
            self.model = f"models/{gene}.clintegrate"

            df = self.load_df(
                f"precomputed_snp_values/cleaned_" + gene + ".csv"
            ).set_index("name")

            df.drop(
                columns = list(
                    filter(lambda x : "Unnamed"  in x, list(df.columns))
                ),
                inplace = True
            )
            self.precomputed_snps = df
            self.initialized = True
        else:
            raise(InvalidGeneException(gene))

    def validate_initialized(self):
        if not self.initialized:
            raise(NotInitializedException())

    def help(self):
        print(help_statement)

    def load_example_data(self):
        self.validate_initialized()
        df = self.load_df(f"example_data/{self.gene}.csv", use_index_0=True)
        return df

    def validate_input_csv(self, df):
        self.validate_initialized()
        required_fields = gene_to_required_fields[self.gene]
        for r in required_fields:
            if r not in df.columns:
                raise(InvalidFormatError(
                        f"Missing required field in csv header: {r}")
                )

        if "sex" in required_fields:
            vals = set(df["sex"])
            for c in vals:
                if c not in {"M", "F"}:
                    raise(InvalidFormatError(f"""
                    Invalid value for sex, allowed valuesare M, F - found {c}
                    """))

        for v in df["PRS"].values:
            try:
                float(v)
            except:
                raise(InvalidFormatError(f"""
                PRS values must be type float, found {v}
                """))

        for v in df["Family History"].values:
            if v not in {0, 1}:
                raise(InvalidFormatError(f"""
                Family History values must ne binary 0/1, found {v}
                """))

    def onset_by_age(df):
        self.validate_initialized()
        df = clean_and_fill_missing_values(df)
        pass

    def clean_input_csv(self, df):
        self.validate_initialized()
        sex_to_int = {"M" : 1, "F": 0}
        if "sex" in df.columns:
            df["sex"] = df["sex"].apply(lambda x : sex_to_int[x])
        return df

    def append_variant_info(self, df):
        self.validate_initialized()
        variant_fields = default_variant_info_dictionaries[self.gene]
        precomputed_variants = set(self.precomputed_snps.index)
        patient_variant_dict = {}
        for index, row in df.iterrows():
            v = row["variant"]
            if type(row["variant"]) == str and row["variant"] != "":
                if v in precomputed_variants:
                    variant_values = self.precomputed_snps.loc[v]
                    patient_variant_dict[index] = variant_values.to_dict()
                else:
                    print("not found")
            else:
                patient_variant_dict[index] = variant_fields

        variant_level_df = pd.DataFrame.from_dict(
            patient_variant_dict,
            orient = "index"
        )

        joined_df = df.join(variant_level_df)
        stream = pkg_resources.resource_stream(__name__, self.model)
        cph = pickle.load(stream)
        predictions = cph.predict_partial_hazard(joined_df)
        predictions.rename(
            columns = {0 : "Partial Hazard Prediction"},
            inplace = True
        )

        predictions_appended = df.join(predictions)
        return predictions_appended

    def validate_reference_alleles(self, variants):
        self.validate_initialized()
        for variant in variants:
            if type(variant) == str and variant != "":
                validate_variant_format(variant, self.gene)



    def generate_risk_predictions(self, df):
        self.validate_initialized()
        self.validate_input_csv(df)
        df = self.clean_input_csv(df)
        self.validate_reference_alleles(df["variant"].values)
        appended_df = self.append_variant_info(df)
        int_to_sex = {1:"M", 0:"F"}
        if "sex" in appended_df.columns:
            appended_df["sex"] = appended_df["sex"].apply(lambda x : int_to_sex[x])
        return appended_df

if __name__ == "__main__":
    ipm = IntegrativePredictiveModel()
    ipm.initialize("APOB")
    data = ipm.load_example_data()
    res = ipm.generate_risk_predictions(data)
    print(res)
