# How to prepare your variant call sets (VCF) for VariantSurvival

### Table of Contents
1. [Preamble](#Preamble)
2. [Merge VCFs across samples](#Merge-VCFs-across-samples)
3. [Annotate genes in the multi-sample VCF file](#Annotate-genes-in-the-multi\-sample-VCF-file)
<br>

## Preamble 

VariantSurvival is a platform to analyze genotype-treatment response with respect to structural variants.
Hence, the foundation for the survival analysis in the VariantSurvival Shiny App are structural variants (SV) and it is expected that an SV callset is present for every individual of the study group.
VariantSurvival is sequencing platform independent and works with callsets of arbitrary variant calling algorithms as long as they comply with the standard VCF file format.
The following sections provide recommendations on how to prepare the input data for the VariantSurvival from multiple SV callsets.

## Merge VCFs across samples

The first requirement to prepare the input data is to merge SV callsets across multiple individuals, e.g. of a study group.
For a clinical study, this must be all individuals of the drug _and_ placebo group.
Merging SV arcoss multiple individuals can be done with a variety of bioinformatics tools (e.g. [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), [Jasmine](https://github.com/mkirsche/Jasmine) or [svimmer](https://github.com/DecodeGenetics/svimmer)).
The merging procedure below uses SURVIVOR and assumes the following generic folder structure in a working directory `work_dir`:
```
work_dir
└── variants
    ├── sample1.vcf
    ├── sample2.vcf
    ├── sample3.vcf
    ├── sample4.vcf
    └── sample5.vcf
```
First, SURVIVOR requires a list of VCF files to merge.
In the folder structure above the command 
```
ls work_dir/variants/*.vcf > work_dir/samples
```
lists all VCF files and stores that list in a new file `samples`.
Second, the SURVIVOR command that actually merges the variant callsets into one multi-sample VCF file could be
```
SURVIVOR merge work_dir/samples 1000 1 1 1 0 30 samples_merged.vcf
```
For more details on the individual parameters please see the official SURVIVOR [Wiki](https://github.com/fritzsedlazeck/SURVIVOR/wiki).
In the next step we will look into how the individual SV within `samples_merged.vcf` can be annotated with gene identifiers in order for VariantSurvival to filter variants of interest.

## Annotate genes in the multi-sample VCF file

_TODO : gff overlap right now? GeneVar soon?_
