# How to prepare your variant call sets (VCF) for VariantSurvival

### Table of Contents
1. [Preamble](#Preamble)
2. [Merge VCFs across samples](#Merge-VCFs-across-samples)
3. [Annotate genes in the multi-sample VCF file](#Annotate-genes-in-the-multi\-sample-VCF-file)
<br>

## Preamble 

VariantSurvival is a platform to analyze genotype-treatment response with respect to structural variants.
Hence, the foundation for the survival analysis in the VariantSurvival Shiny App are structural variants (SV) and it is expected that an SV callset is present for every individual of a study group.
VariantSurvival is sequencing platform independent and works with callsets of arbitrary variant calling algorithms as long as they comply with the standard VCF file format.
The following sections provide recommendations on how to prepare the input data for the VariantSurvival from multiple SV callsets.

## Merge variant callsets across individuals

### TL;DR:
Use a tool of your choice to merge SV callsets across all individuals of your study groups. Some recommendations listed below:
- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
- [Jasmine](https://github.com/mkirsche/Jasmine)
- [svimmer](https://github.com/DecodeGenetics/svimmer)

### Step-by-step:
The first step towards the input file format of VariantSurvival is to merge SV callsets across multiple individuals, e.g. of a study group.
For a clinical study, this must be all individuals of the study, i.e. of the drug _and_ placebo groups.
Merging SV arcoss multiple individuals can be done with a variety of bioinformatics tools (see recommendations listed above).
Here, we demonstrate the merging procedure using SURVIVOR and assume the following generic folder structure:
```
project
└── variants
    ├── sample1.vcf
    ├── sample2.vcf
    ├── sample3.vcf
    ├── sample4.vcf
    └── sample5.vcf
```
First, SURVIVOR requires a list of VCF files to merge.
From your shell you can list all your VCF files and write them into a text file called _samples_ via
```
ls variants/*.vcf > samples
```
Next, we use SURVIVOR to actually merge the variant callsets (now listed in the 'samples' file) into one multi-sample VCF file via
```
SURVIVOR merge samples 1000 1 1 1 0 30 samples_merged.vcf
```
For more details on the individual parameters of the `SURVIVOR merge` command please see the official SURVIVOR [Wiki](https://github.com/fritzsedlazeck/SURVIVOR/wiki).
In the next step we will look into how each individual SV within `samples_merged.vcf` can be annotated with gene identifiers in order for VariantSurvival to filter variants of interest.

## Annotate genes in the multi-sample VCF file

### TL;DR:
Use a tool of your choice to annotate the variants of your study group with gene identifiers. Some recommendations listed below:
- [bcftools](https://samtools.github.io/bcftools/)
- [vcfanno](https://github.com/brentp/vcfanno)
- [vcf-annotator](https://github.com/rpetit3/vcf-annotator)

### Step-by-step
The first step towards the input file format of VariantSurvival is to annotate the merged SV callset of your study group with gene identifiers.
A list of gene identifiers together with their genomic coordinates can be retrieved from arbitrary databases or resources.
The format of that list depends on the annotation tool of choice.
We will demonstrate the annotation procedure using bcftools, hence we require the list in the standard [BED file](http://www.ensembl.org/info/website/upload/bed.html) format.
Given the minimum information required in the leading four columns of a BED file `genes.bed.gz`, like
| CHROM | FROM | TO | GENE |
| --- | --- | --- | --- |
| chr1 | 11200 | 11500 | TERC |
| chr1 | 465076 | 465431 | SOX |
| chr4 | 5333101 | 5333440 | HOXD |
| ... | ... | ... | ... |

we can use bcftools to annotate the variant records of our merged VCF file via
```
bgzip -c samples_merged.vcf > samples_merged.vcf.gz
tabix -p vcf samples_merged.vcf.gz
bcftools annotate \
  -a genes.bed.gz \
  -c CHROM,FROM,TO,GENE \
  -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  samples_merged.vcf.gz
```
In more detail, the first two commands prepare (compress and index) the VCF file resulting from the merge step.
The `bcftools annotate` command adds the gene identifiers to the variant records' INFO field and modifies the VCF file's header section accordingly.
