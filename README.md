# VariantSurvival: A tool to identify genotype-treatment response
<img src="https://user-images.githubusercontent.com/41301333/195215088-8404f200-8297-4322-a30f-c84f526aa620.png" width="200" height="200" align="right">

### Table of Contents
1. [About VariantSurvival](#about-variantsurvival)
2. [The developers team](#the-developers-team)
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

The second strictly required input file...

VariantSurvival run example:
```
library(VariantSurvival)
VariantSurvival::VariantSurvival(vcffile="my_variant_file.vcf", metadatafile= "my_metadata_file.xlsx")
```


## The developers team

* Ahmad Al Khleifat [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/AhmadAlKhleifat.svg?style=social&label=Follow%20%40AhmadAlKhleifat)](https://twitter.com/AhmadAlKhleifat)
* Thomas Krannich [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/krannich479.svg?style=social&label=Follow%20%40krannich479)](https://twitter.com/krannich479)
* Hiba Ben Aribi [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/Hiba_BenAribi.svg?style=social&label=Follow%20%40Hiba_BenAribi)](https://twitter.com/Hiba_BenAribi)
* Marina Herrera Sarrias
* Moustafa Shokrof [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/mostafashokrof2.svg?style=social&label=Follow%20%40mostafashokrof2)](https://twitter.com/mostafashokrof2)


