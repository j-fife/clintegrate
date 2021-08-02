included_genes = [
    "BRCA1",
    "BRCA2",
    "LDLR",
    "PCSK9",
    "APOB"
]

GENE_TO_START_38 = {
    "APOB" : 21001429,
    "BRCA1" : 43044295,
    "BRCA2" : 32315086,
    "LDLR" : 11089432,
    "PCSK9" : 55039447
}

GENE_TO_END_38 = {
    "APOB" : 21044073,
    "BRCA1" : 43170245,
    "BRCA2" : 32400268,
    "LDLR" : 11133820,
    "PCSK9" : 55064852
}

GENE_TO_CHR = {
    "APOB" : 2,
    "BRCA1" : 17,
    "BRCA2" : 13,
    "LDLR" : 19,
    "PCSK9" : 1
}

default_variant_info_dictionaries = {
    "LDLR" : {
        "CADD"  : 0,
        "GERP"  : 0,
        "phyloP" : 0,
        "Synonymous" : 0,
        "Deleterious"  : 0,
        "Missense" : 0,
        "allele_frequency" : -2,
        "Region 1" : 0,
        "Region 2" : 0,
        "Region 3" : 0,
        "Region 4" : 0,
        "Region 5" : 0
    },
    "APOB" : {
        "CADD"  : 0,
        "GERP"  : 0,
        "phyloP" : 0,
        "Synonymous" : 0,
        "Deleterious"  : 0,
        "Missense" : 0,
        "allele_frequency" : -2,
        "Region 1" : 0,
        "Region 2" : 0,
        "Region 3" : 0,
        "Region 4" : 0,
        "Region 5" : 0,
        "Region 6" : 0,
        "Region 7" : 0,
        "Region 8" : 0,
        "Region 9" : 0,
        "Region 10" : 0,
        "Region 11" : 0,
        "Region 12" : 0,
        "Region 13" : 0
    },
    "BRCA1" :{
        "CADD"  : 0,
        "GERP"  : 0,
        "phyloP" : 0,
        "Synonymous" : 0,
        "Deleterious"  : 0,
        "Missense" : 0,
        "allele_frequency" : -2,
        "Region 1" : 0,
        "Region 2" : 0,
        "Region 3" : 0,
        "Region 4" : 0,
        "Region 5" : 0
    },
    "BRCA2" :{
        "CADD"  : 0,
        "GERP"  : 0,
        "phyloP" : 0,
        "Synonymous" : 0,
        "Deleterious"  : 0,
        "Missense" : 0,
        "allele_frequency" : -2,
        "Region 1" : 0,
        "Region 2" : 0,
        "Region 3" : 0,
        "Region 4" : 0,
        "Region 5" : 0,
        "Region 6" : 0,
        "Region 7" : 0,
        "Region 8" : 0,
        "Region 9" : 0,
        "Region 10" : 0,
    },
    "PCSK9" : {
        "CADD"  : 0,
        "GERP"  : 0,
        "phyloP" : 0,
        "Synonymous" : 0,
        "Deleterious"  : 0,
        "Missense" : 0,
        "allele_frequency" : -2,
        "Region 1" : 0,
        "Region 2" : 0,
        "Region 3" : 0,
        "Region 4" : 0,
        "Region 5" : 0,
    }


}

gene_to_region_breaks = {
    "APOB" : [0, 819, 1279, 2604, 4089, 5238, 6835, 7939, 8392, 9491, 10187, 10808, 11789, 13220],
    "BRCA1": [0, 928, 2215, 3514, 4689, 5579],
    "BRCA2": [0, 979, 2151, 3328, 4444, 5453, 6395, 7375, 8219, 8915, 10234],
    "LDLR":  [0, 465, 895, 1394, 1916, 2548],
    "PCSK9": [0, 426, 658, 1274, 1503]
}


GENE_TO_CONDITION = {
    "APOB" : "Coronary Artery Disease",
    "BRCA1" : "Breast Cancer",
    "BRCA2" : "Breast Cancer",
    "LDLR" : "Coronary Artery Disease",
    "PCSK9" : "Coronary Artery Disease"
}

gene_to_required_fields = {
    "BRCA1" : ["PRS", "Family History", "variant"],
    "BRCA2" : ["PRS", "Family History", "variant"],
    "LDLR" : ["PRS", "Family History", "sex", "variant"],
    "APOB" : ["PRS", "Family History", "sex", "variant"],
    "PCSK9" : ["PRS", "Family History", "sex", "variant"]
}

gene_string = ", ".join(included_genes)



gene_to_assumptions = {

}

gene_to_info = {

}

gene_to_example_data = {

}

version = "0.0.1"

help_statement = f"""
____________________________________________________________________
      _ _       _                       _
  ___| (_)_ __ | |_ ___  __ _ _ __ __ _| |_ ___
 / __| | | '_ \| __/ _ \/ _` | '__/ _` | __/ _ |
| (__| | | | | | ||  __/ (_| | | | (_| | ||  __/
 \___|_|_|_| |_|\__\___|\__, |_|  \__,_|\__\___|
                        |___/

Version {version}
____________________________________________________________________

Welcome to clintegrate!

To get started initialize a model, and load the example data:

TODO:

Not all models use the same data, so be sure to check each time you
use a new gene's predictive model.

See further documentation here:
TODO

If you used Clintegrate in your projects please consider citing us:
TODO

"""
