from constants import (
    default_variant_info_dictionaries,
    gene_to_region_breaks,
    GENE_TO_END_38,
    GENE_TO_START_38,
    GENE_TO_CHR
)
from Exceptions import (
    VariantNotInGeneException,
    VariantFormatException,
    ReferenceSequenceException
)
import requests
import json
import numpy as np
import pandas as pd

EXT = '/?CADD=1&canonical=1&dbNSFP=phyloP100way_vertebrate,GERP%2B%2B_RS&content-type=application/json'
VEP38_URL = 'https://rest.ensembl.org/vep/human/hgvs/'

def determine_most_severe_consequence(terms_string):
    delet_terms = ["frameshift","stop_lost","splice_donor","splice_acceptor","stop_gain", "nonsense"]
    for term in delet_terms:
        if term in terms_string:
            return "Delet-" + term
    if "missense" in terms_string:
        return "Missense"
    if "synon" in terms_string:
        return "Synonymous"

    if "inframe_del" in terms_string or "inframe_ins" in terms_string:
        return "Inframe Indel"

    else:
        return "Non-Coding"

def place_position(breaks, pos):
    i = 0
    while pos > breaks[i] and i < len(breaks) - 1:
        i += 1
    return i

def create_usable_dict(d, gene):
    ## Need to fix regions
    to_return = default_variant_info_dictionaries[gene]
    region_breaks = gene_to_region_breaks[gene]
    n_breaks = len(region_breaks)-1
    for i in range(1, n_breaks + 1):
        to_return["Region " + str(i)] = 0
    if "Delet" in d["consequence"]:
        to_return["Deleterious"] = 1
        to_return["CADD"] = d["CADD"]
        to_return["GERP"] = d["GERP"]
        to_return["phyloP"] = d["phyloP"]
        to_return["allele_frequency"] = d["af"]

    elif d["consequence"] == "Missense":
        to_return["Missnese"] = 1
        to_return["CADD"] = d["CADD"]
        to_return["GERP"] = d["GERP"]
        to_return["phyloP"] = d["phyloP"]
        region = place_position(region_breaks, d["cds"])
        to_return["Region " + str(region)] = 1
        to_return["allele_frequency"] = d["af"]

    elif d["consequence"] == "Synonymous":
        to_return["Synonymous"] = 1
        to_return["CADD"] = d["CADD"]
        to_return["GERP"] = d["GERP"]
        to_return["phyloP"] = d["phyloP"]
        to_return["allele_frequency"] = d["af"]

    return to_return


def parse_vep_json(j, gene):
    canonical_info = None
    for val in j[0]["transcript_consequences"]:
        if "canonical" in val.keys():
            canonical_info = val
    af = -6
    if "colocated_variants" in j[0].keys():
        if len(j[0]["colocated_variants"]) > 0:
             coloc_dict = j[0]["colocated_variants"][0]
             if type(coloc_dict) == dict and "frequencies" in coloc_dict.keys():
                 freq_dict = coloc_dict["frequencies"]
                 if len(freq_dict) > 0:
                     vals = freq_dict[list(freq_dict.keys())[0]].values()
                     max_af = np.max(list(vals))
                     if max_af > 0:
                         af = np.log10(max_af)
    cds = canonical_info.get("cds_start")
    CADD = canonical_info.get("cadd_raw", 0)
    gerp = canonical_info.get("gerp++_rs", 0)
    phylop = canonical_info.get("phylop100way_vertebrate", 0)
    cosequence_terms = canonical_info.get("consequence_terms", [])

    if len(cosequence_terms) > 0 :
        consequence = determine_most_severe_consequence(",".join(cosequence_terms))
    else:
        consequence = "Non-Coding"

    d = {
        "consequence" :  consequence,
        "cds" : cds,
        "CADD" : CADD,
        "GERP" : gerp,
        "phyloP" : phylop,
        "af" : af
    }

    return create_usable_dict(d, gene)


def deletion_get_variant_data(chrom, start, length, gene):
    variant = str(chrom) + ":g." + str(start) + "_" + str(start + length) + "del"
    url = VEP38_URL
    res = requests.get(url + variant + EXT)
    print(res.json())
    return parse_vep_json(res.json(), gene)


def insertion_get_variant_data(chrom, start, inserted_sequence, gene):
    # print(start + len(inserted_sequence))
    variant = str(chrom) + ":g." + str(start) + "_" + str(start + 1) +  "ins" + inserted_sequence
    url = VEP38_URL
    res = requests.get(url + variant + EXT)
    return parse_vep_json(res.json(), gene)

def SNP_get_variant_data(chr, start, ref, alt, gene):
    variant = str(chr) + ':g.' + str(start) + str(ref) + '>' + str(alt)
    url = 'https://rest.ensembl.org/vep/human/hgvs/'
    res = requests.get(url + variant + EXT)
    # canonical_transcript = [t for t in data.get('transcript_consequences', []) if 'canonical' in t and t['gene_symbol']==gene]
    return parse_vep_json(res.json(), gene)
    #
    # except:
    #     print("threw exception when getting json from vep" , url + variant + EXT)


def validate_variant_format(variant, gene):
    start = GENE_TO_START_38[gene]
    end = GENE_TO_END_38[gene]
    chrom = GENE_TO_CHR[gene]

    with open(f'./data/sequences/{gene}_38.json') as f:
        gene_json = json.load(f)
    f.close()

    approved_characters = set(["A", "T", "C", "G"])
    try:
        pieces = variant.split("-")
    except:
        raise(VariantFormatException(variant))


    if len(pieces) != 4:
        raise(VariantFormatException(variant))
    if pieces[0] != str(chrom):
        raise(VariantNotInGeneException(variant, chrom, start, end, gene))

    try:
        pos = int(pieces[1])
        if pos >= start and pos <= end:
            for char in pieces[2]:
                if char not in approved_characters:
                    raise(VariantFormatException(variant))
            for char in pieces[3]:
                if char not in approved_characters:
                    raise(VariantFormatException(variant))
        else:
            raise(VariantNotInGeneException(variant, chrom, start, end, gene))
    except:
        raise(VariantFormatException(variant))


    if len(pieces[2]) == 1 and len(pieces[3]) == 1:
        sequence = gene_json["sequence"]
        difference = pos - gene_json["start"]
        ref_allele = sequence[difference]
        if ref_allele != pieces[2]:
            raise(ReferenceSequenceException(ref_allele, pieces[2], variant))

def clean_and_fill_missing_values(df, gene):
    gene_start = None
    gene_end = None
    patient_to_variant_info_dicts = {}
    default_dict = default_variant_info_dictionaries[gene]
    precomputed_variants_file = pd.read_csv(f"./precomputed_snp_values/{gene}.csv")
    for index, row in df.iterrows():
        variant = row["variant"]
        if variant != "" and variant is not None:
            validate_variant_format(variant)
            pieces = variant.split("-")
            position = int(pieces[1])
            chrom = int(piees[0])
            if len(pieces[2]) == 1 and len(pieces[3]) == 1:
                data = precomputed_variants_file.loc[variant].to_dict()
            elif len(pieces[2]) > len(pieces[3]):
                data = deletion_get_variant_data(
                    chrom,
                    position + len(pieces[3]) - len(pieces[2]),
                    len(pieces[3]) - len(pieces[2])
                )
            elif len(pieces[3]) > len(pieces[2]):
                if len(pieces[2]) != 1 or pieces[3][0] != pieces[2]:
                    raise(VariantFormatException(variant))
                data = insertion_get_variant_data(chrom, position + 1, pieces[3][1:])
            else:
                data = default_dict
        else:
            data = default_dict
        patient_to_variant_info_dicts[index] = data

    variant_info_df = pd.DataFrame(patient_to_variant_info_dicts)
    df = df.join(variant_info_df)
    return df
    ## need to clean df


if __name__ == "__main__":
    #validate_variant_format("2-21001433-T-G", "APOB")
    i = deletion_get_variant_data(2, 21001735, 1, "APOB")
    print(i)
    k = SNP_get_variant_data(2, 21001735, "G", "A", "APOB")
    print(k)
