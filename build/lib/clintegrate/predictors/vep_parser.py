import pandas as pd
import copy
import re
import numpy as np
from clintegrate.predictors.constants import default_variant_info_dictionaries, gene_to_region_breaks

def annotate_consequences(consequence_list):
    consequences_used = []
    is_coding = []
    ## getting varaint type and functional consequence
    for c in consequence_list:
        if "splice_acceptor" in c or "frameshift_variant" in c or "splice_donor_variant" in c or "start_lost" in c or "stop_gained" in c:
            consequences_used.append("Deleterious")
            is_coding.append(True)
        elif "missense" in c:
            consequences_used.append("Missense")
            is_coding.append(True)
        elif "inframe_insertion" in c or "inframe_deletion" in c:
            consequences_used.append("Inframe Indel")
            is_coding.append(True)
        elif "synon" in c:
            consequences_used.append("Synonymous")
            is_coding.append(False)
        else:
            consequences_used.append("Intronic")
            is_coding.append(False)
    return consequences_used, is_coding


def create_input_df(variant_list):
    header = ["#CHROM", "POS","ID","REF", "ALT", "QUAL", "FILTER", "INFO"]
    rows = []
    for v in variant_list:
        pieces = v.split("-")
        chrom = int(pieces[0])
        pos = int(pieces[1])
        ref = pieces[2]
        alt = pieces[3]
        rows.append([chrom, pos, ".", ref, alt, ".", "PASS", "."])
    vcf_df = pd.DataFrame(rows, columns = header)
    vcf_df.to_csv(f"./data/vcf_files/{gene}_SNPs.vcf", sep = "\t", index = False)


def parse_annotated_vep_file(
    input_file
):
    rows_to_skip = 0
    headers = []
    seen_VEP_fields = set()
    with open(input_file, "r") as f:
        found = False
        line = True
        while not found and line:
            line = f.readline()
            if "#CHROM" in line and "REF" in line:
                found = True
                break
            rows_to_skip += 1
            headers.append(line)
    f.close()
    df = pd.read_csv(input_file, delimiter = "\t", skiprows = rows_to_skip)
    info_sections = headers[2].split("Ensembl VEP. Format: ")[1].replace(">\n", "").split("|")[1:]
    columns = ["Name", "Chrom", "Pos", "Ref", "Alt"] + info_sections
    rows = []

    ## 399 fields are returned by VEP with dbNSFP plugin. run using "ALL" as
    ## field filter instead of
    n = 399
    for index, row in df.iterrows():
        chrom, pos, ref, alt = row["#CHROM"], int(row["POS"]), row["REF"], row["ALT"]
        if chrom not in ["X", "Y"]:
            chrom = int(chrom)
        name = str(chrom) + "-" + str(pos) + "-" + ref + "-" + alt
        pieces =  list(map(lambda x: None if x == "" else x, row["INFO"].split("|")[1:]))
        # assert len(pieces) % n == 0
        transcripts = [pieces[i:i + n] for i in range(0, len(pieces), n)]
        canonical_fields = None
        canonical_found = 0
        for c in transcripts:
            d = {}
            for field, value in zip(info_sections, c):
                d[field] = value

            if d["CANONICAL"] and "YES" in d["CANONICAL"]:
                canonical_fields = c
                canonical_found = 1
                break

        if canonical_fields is None:
            canonical_fields = transcripts[0]

        assert len(pieces) % n == 0
        rows.append([name, chrom, pos, ref, alt] + canonical_fields + [canonical_found])

    variant_df = pd.DataFrame(rows, columns = columns + ["canonical_found"])
    consequences_used, is_coding = annotate_consequences(variant_df["Consequence"].values)
    variant_df["consequence_use"] = consequences_used
    variant_df["is_coding"] = is_coding

    # getting  continuous values
    ## AF will eventually be overwritten from REST VEP calls
    variant_df["af_use"] = list(map(lambda x : x if x is None else float(x), variant_df["gnomAD_genomes_POPMAX_AF"]))
    variant_df["CADD_use"] = list(map(lambda x : x if x is None else float(x), variant_df["CADD_raw"]))
    variant_df["GERP_use"] = list(map(lambda x : x if x is None else float(x), variant_df["GERP++_RS"]))
    variant_df["phyloP_use"] = list(map(lambda x : x if x is None else float(x), variant_df["phyloP100way_vertebrate"]))
    variant_df["AA_ref_use"] = list(map(lambda x : "Stop" if x == "X" else x, variant_df["aaref"]))
    variant_df["AA_alt_use"] = list(map(lambda x : "Stop" if x == "X" else x, variant_df["aaalt"]))
    variant_df["AA_pos_use"] = list(map(lambda x : x if x is None else int(list(filter(lambda x : x != "", re.split(r'\D+', x)))[0]), variant_df["aapos"]))
    variant_df["coding_position_use"] = list(map(lambda x : x if x is None else int(list(filter(lambda x : x != "", re.split(r'\D+',x)))[0]), variant_df["CDS_position"]))
    ## need CpG context and transition + transversion filters

    columns_keeping = ['Name', 'Chrom', 'Pos', 'Ref','Alt', "SYMBOL" ,'Consequence',"VEP_canonical",
                    "consequence_use","is_coding", "af_use",  "CADD_use", "GERP_use", "phyloP_use",
                     "AA_ref_use", "AA_alt_use", "AA_pos_use", "coding_position_use",
                     ]
    all_columns = list(variant_df.columns)
    for c in columns_keeping:
        all_columns.remove(c)

    variant_df = variant_df.drop(columns = all_columns)
    variant_df.rename(columns = {
        "Consequence" : "vep_consequence",
        "consequence_use" : "consequence",
        "af_use" : "allele_frequency",
        "SYMBOL" : "Symbol",
        "CADD_use" : "CADD",
        "GERP_use" : "GERP",
        "phyloP_use" : "phyloP",
        "AA_ref_use" : "AA_ref",
        "AA_alt_use" : "AA_alt",
        "AA_pos_use" : "AA_pos",
        "coding_position_use" : "coding_position"
    }, inplace = True)

    variant_df = variant_df.set_index("Name")
    return variant_df


def place_position(breaks, pos):
    i = 0
    while pos > breaks[i] and i < len(breaks) - 1:
        i += 1
    return i

def create_precomputed_SNP_df(input_annotated_vcf, gene):
    variant_df = parse_annotated_vep_file(input_annotated_vcf)
    regions = gene_to_region_breaks[gene]
    default_dict = default_variant_info_dictionaries[gene]
    rows = []
    for index, row in variant_df.iterrows():
        name = row.name
        d = copy.deepcopy(default_dict)
        d["name"] = name
        if row["consequence"] == "Deleterious":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af
            d["Deleterious"] = 1

        elif row["consequence"] == "Missense":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            d["Missense"] = 1
            region = place_position(regions, row["coding_position"])
            d[f"Region {region}"] = 1
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af
        elif row["consequence"] == "Synonymous":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            d["Synonymous"] = 1
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af

        rows.append(d)

    to_return = pd.DataFrame(rows)
    return to_return
