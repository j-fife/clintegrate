          _ _       _                       _      
      ___| (_)_ __ | |_ ___  __ _ _ __ __ _| |_ ___
     / __| | | '_ \| __/ _ \/ _` | '__/ _` | __/ _ |
    | (__| | | | | | ||  __/ (_| | | | (_| | ||  __/
     \___|_|_|_| |_|\__\___|\__, |_|  \__,_|\__\___|
                            |___/

Version 0.0.1
____________________________________________________________________

![PyPI Latest Release](https://img.shields.io/pypi/v/pandas.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)


[Online Web Tool](https://clintegrate.herokuapp.com/)

# Getting Started


The following disease models are available with varaint integration in genes associated with each condition:

|        Condition        	|       Genes       	|
|:-----------------------:	|:-----------------:	|
| Breast Cancer           	| BRCA1, BRCA2      	|
| Colorectal Cancer 	| MLH1, MSH2, MSH6, PMS2 	|
| Coronary Artery Disease 	| APOB, LDLR, PCSK9 	|


#### Initalizing a model

```python
from clintegrate.predictors import IntegrativePredictiveModel as ipm

# Making predictions of coronary artery disease risk using LDLR variants
ipm.initialize("LDLR")
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

