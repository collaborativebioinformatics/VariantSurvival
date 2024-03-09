# VariantSurvival: A tool to identify genotype-treatment response
<img src="https://user-images.githubusercontent.com/41301333/195215088-8404f200-8297-4322-a30f-c84f526aa620.png" width="200" height="200" align="right">

### Table of Content
1. [About VariantSurvival](#about-variantsurvival) <br>
1.1 [Installation](#installation) <br>
1.2 [Running the demo](#running-the-demo) <br>
1.3 [Running your data analysis](#running-your-data-analysis) <br>
2. [Citing](#citing)
3. [The developers team](#the-developers-team)
<br>

## About VariantSurvival
We present VariantSurvival, a lightweight application to visualize genotype-treatment response.
The R application launches a dashboard to browse variant abundances and survival statistics for neurological diseases and corresponding target genes.
VariantSurvival requires as input a set of annotated structural variants in the standard Variant Call Format (VCF) and a file of tabular metadata about the clinical trial or cohorts.
<br>

### Installation
The application is installed via the R programming language using the following commands.
First, the _devtools_ R package is required to download VariantSurvival from its codebase.
If you don't have devtools installed in your active R environment use
```R
install.packages('devtools') #install devtools package
library(devtools) #activate devtools package
```

Now VariantSurvival itself can be installed (downloaded and activated) via
```R
devtools::install_github("collaborativebioinformatics/VariantSurvival/VariantSurvival_package")
library(VariantSurvival)
```

Note that all dependencies will be installed automatically when installing VariantSurvival.
After a successful installation the currently installed software version can be retrieved via
```R
packageVersion("VariantSurvival")
```

### Running the demo
After following the installation instructions the VariantSurvival dashboard can be launched with demo data via
```R
VariantSurvival::VariantSurvival(demo=TRUE)
```

### Running your data analysis
To launch the dashboard with your own clinical data two input files have to be provided.
One required input is a VCF file with a set of gene-annotated structural variants across the entire study cohort.
We have compiled a dedicated manual page on [how to prepare your VCF file](https://github.com/collaborativebioinformatics/VariantSurvival/blob/main/docs/preprocessing/prepareVCF.md) for VariantSurvival with instructions and recommendations how to merge and annotate VCF files. <br>

The second strictly required input file is a tabular metadata Excel sheet.
The minimal set of required columns for the quantitavite and survival analysis is
| Column name    | Type         | Encoding                                  |
| :------------- | :----------- | :---------------------------------------- |
| Patient IDs    | [string]     | Matching the sample names of the VCF file |
| Time factor    | [continuous] | Time since onset in days or years         |
| Trial group    | [binary]     | Treatment / control                       |
| Survival state | [binary]     | Deceased / alive                          |

Given a VCF file `my_variant_file.vcf` and an Excel sheet `my_metadata_file.xlsx` the dashboard is launched passing both files as function paramters.
```R
VariantSurvival::VariantSurvival(vcffile="my_variant_file.vcf", metadatafile= "my_metadata_file.xlsx")
```

## Citing
If you use VariantSurvival, please cite this [journal publication](https://www.frontiersin.org/articles/10.3389/fbinf.2023.1277923/full)

- _Krannich T, Sarrias MH, Ben Aribi H, Shokrof M, Iacoangeli A, Al-Chalabi A, Sedlazeck FJ, Busby B and Al Khleifat A (2023) VariantSurvival: a tool to identify genotype–treatment response. Front. Bioinform. 3:1277923. doi: 10.3389/fbinf.2023.1277923_

## The developers team
* Ahmad Al Khleifat [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/AhmadAlKhleifat.svg?style=social&label=Follow%20%40AhmadAlKhleifat)](https://twitter.com/AhmadAlKhleifat)
* Thomas Krannich [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/krannich479.svg?style=social&label=Follow%20%40krannich479)](https://twitter.com/krannich479)
* Hiba Ben Aribi [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/Hiba_BenAribi.svg?style=social&label=Follow%20%40Hiba_BenAribi)](https://twitter.com/Hiba_BenAribi)
* Marina Herrera Sarrias
* Moustafa Shokrof [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/mostafashokrof2.svg?style=social&label=Follow%20%40mostafashokrof2)](https://twitter.com/mostafashokrof2)

## Help 
For questions about VariantSurvival and bug reports please refer to the GitHub issues section of this repository.
