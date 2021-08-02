          _ _       _                       _      
      ___| (_)_ __ | |_ ___  __ _ _ __ __ _| |_ ___
     / __| | | '_ \| __/ _ \/ _` | '__/ _` | __/ _ |
    | (__| | | | | | ||  __/ (_| | | | (_| | ||  __/
     \___|_|_|_| |_|\__\___|\__, |_|  \__,_|\__\___|
                            |___/

Version Beta 2.1.2
____________________________________________________________________

[![Generic badge](https://img.shields.io/badge/Creator-Christopher_A._Cassa_Lab-maroon.svg)](http://genetics.bwh.harvard.edu/wiki/cassa/)
[![Generic badge](https://img.shields.io/badge/Maintainer-James_Fife-maroon.svg)](https://github.com/j-fife/)

![PyPI Latest Release](https://img.shields.io/pypi/v/pandas.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)
[![Open Source? Yes!](https://badgen.net/badge/Open%20Source%20%3F/Yes%21/blue?icon=github)](https://github.com/Naereen/badges/)
[![DOI:___](https://zenodo.org/badge/DOI/_______.svg)](https://doi.org/___)

[Online Web Tool](https://clintegrate.herokuapp.com/)


# About
Read our full manuscript [Here](./#) (preprint)

Using exome seqencing and data from over 200,000 individuals from the UK Biobank, Clintegrate is a intergrative
risk prediction framework designed to address issues around variant interpretation as it pertains to current
practices in personalized risk assessment.

From our manuscript:
>Mapping germline variants to personalized clinical risk is a major goal in precision medicine.[1] While clinical
>diagnostic testing has advanced dramatically, the interpretation of monogenic variants remains challenging,
>and is generally done at the variant level.[2] To date, such genetic testing has largely been conducted in
>the presence of a phenotypic indication, where the prior probability of detecting a causal variant is high.[3]
>However, subsequent population screening efforts have identified substantial incomplete penetrance or reduced
>expressivity in these previously identified monogenic disease variants[4–6] and assessed how the risk attributable
>to these variants can be modified by clinical and polygenic risk factors.[7–9]

We present a framework that integrates patient-level information
(Polygenic Risk Scores, Family History, Sex, etc.), monogenic variant-level features (CADD, Allele Frequency, GERP,
CpG Context, phyloP, etc.) and protein regional information to create an overall measurement of risk. With growing
availability of population-level sequencing efforts and diagnostic data, we aim to continually expand this model
to include additioanl phenotypes and risk factors (genetic and non-genetic).   

# Getting Started


The following disease models are available with varaint integration in genes associated with each condition:

|        Condition        	|       Genes       	|
|:-----------------------:	|:-----------------:	|
| Breast Cancer           	| BRCA1, BRCA2      	|
| Colorectal Cancer 	| MLH1, MSH2, MSH6, PMS2 	|
| Coronary Artery Disease 	| APOB, LDLR, PCSK9 	|


&nbsp;
&nbsp;
&nbsp;

## Initalizing a model

```python
from clintegrate.predictors.risk import IntegrativePredictiveModel

# Making predictions of coronary artery disease risk using APOB variants
ipm = IntegrativePredictiveModel()
ipm.initialize("APOB")
```

Disease models have different required fields for accurate risk assessments, which
we're actively updating. It's always a good idea to look at the example data:

```python
ipm.load_example_data()
```
| id      | sex   |   PRS |   Family History | variant        |
|:--------|:------|------:|-----------------:|:---------------|
| person1 | M     | -0.54 |                0 |                |
| person2 | F     |  2.51 |                1 | 2-21001432-G-A |
| person3 | F     |  0    |                1 | 2-21001769-G-T |
| person4 | F     |  1.3  |                0 |                |


Variants are always optional for making risk assessments, however the field and
the remaining fields are required.  



&nbsp;
&nbsp;
&nbsp;

## Making predictions
To predict [partial hazard](https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html#cox-s-proportional-hazard-model) values

```python
example_data = ipm.load_example_data()
generate_risk_predictions(example_data)
```

| id      |   sex |   PRS |   Family History | variant        |   Partial Hazard Prediction |
|:--------|------:|------:|-----------------:|:---------------|----------------------------:|
| person1 |     M | -0.54 |                0 |                |                    1.37421  |
| person2 |     F |  2.51 |                1 | 2-21001432-G-A |                    2.06113  |
| person3 |     F |  0    |                1 | 2-21001769-G-T |                    0.230054 |
| person4 |     F |  1.3  |                0 |                |                    0.645135 |


&nbsp;
&nbsp;
&nbsp;

### References


[1] Green, E. D. et al. Strategic vision for improving human health at The Forefront of Genomics. Nature 586, 683–692 (2020).

[2]	Richards, S. et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet. Med. 17, 405–423 (2015).

[3]	Rehm, H. L. et al. ClinGen — The Clinical Genome Resource. N. Engl. J. Med. 372, 2235–2242 (2015).

[4]	Goodrich, J. K. et al. Determinants of penetrance and variable expressivity in monogenic metabolic conditions across 77,184 exomes. Nat. Commun. 12, 3505 (2021).

[5]	Cassa, C. A., Tong, M. Y. & Jordan, D. M. Large numbers of genetic variants considered to be pathogenic are common in asymptomatic individuals. Hum Mutat (2013) doi:10.1002/humu.22375.

[6]	Forrest, I. S. et al. Ancestrally and Temporally Diverse Analysis of Penetrance of Clinical Variants in 72,434 Individuals. http://medrxiv.org/lookup/doi/10.1101/2021.03.11.21253430 (2021) doi:10.1101/2021.03.11.21253430.

[7]	Fahed, A. C. et al. Polygenic background modifies penetrance of monogenic variants conferring risk for coronary artery disease, breast cancer, or colorectal cancer. http://medrxiv.org/lookup/doi/10.1101/19013086 (2019) doi:10.1101/19013086.

[8]	Friebel, T. M., Domchek, S. M. & Rebbeck, T. R. Modifiers of Cancer Risk in BRCA1 and BRCA2 Mutation Carriers: A Systematic Review and Meta-Analysis. JNCI J. Natl. Cancer Inst. 106, (2014).

[9]	Saadatagah, S. et al. Genetic basis of hypercholesterolemia in adults. Npj Genomic Med. 6, 1–7 (2021).
