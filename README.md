# VariantSurvival : A tool to identify genotype-treatment response
<img src="https://user-images.githubusercontent.com/41301333/195215088-8404f200-8297-4322-a30f-c84f526aa620.png" width="300" height="300">


##  Abstract

For a number of neurological diseases, such as Alzheimer's disease, Parkinson's disease and many others, certain genes are known to be involved in the disease mechanism. A common question is whether a structural variant in any such gene may be related to drug response in clinical trials, and how this relationship can contribute to the lifecycle of drug development. To this end, we introduce VariantSurvival, a tool that identifies changes in survival relative to structural variants within target genge. VarantSurvival matches annotated structural variants with clinically relevant genes of neurological disease. A cox regression model determines the change in survival between placebo and clinical trial group with respect to the number of structural variants in the drug target genes. We showcase the functionality of our approach on the example of Amyotrophic lateral sclerosis (ALS) and the SETX gene. VariantSurvivor has a user-friendly and lightweight graphical user interface built on the shiny webapp package.

##  Goal and motivation
The genomics data-driven identification of gene variants has been utilised for predicting survival and clinical outcome such as clinical management and treatments. Many studies have investigated the relationship between gene and survival using large cohorts however, integrating genomics with clinical trials outcome remain underutilized. Here we have developed VariantSurvival a bioinformatics tool evaluates the associations between genomic variant and survival. VariantSurvival is web tool to perform survival analysis on next generation sequencing data (NGS) and a variety of other phenotypic inputs such as age, sex, and other clinical information. VariantSurvival contains multiple genetic information from neurological and psychiatric conditions which will be linked with a list of known genes in the disease of interest. Users can upload NGS data and select gene list to determine the effect on survival outcomes. The tool provides results including box plots of low and high-risk groups, Cox proportional plots in relation to treatment profile. 


##  Introduction

The challenges in developing novel therapeutics for neurodegenerative diseases (ND) result from the paucity of novel, valid targets. This in turn results from etiological heterogeneity, the complex and often polygenic nature of genetic risk. Despite the increase in funding for drug discovery, only 10% of new drug candidates in early stage clinical trials are eventually approved. Recent study has found that drug targets with genetic support were twice as likely to be approved. The discovery of rare and common genetic variants associated with risk for neurological and neuropsychiatric illness provides the opportunity to restart hypothesis-led clinical trials data analysis. High-risk mutations in single genes that identify specific targets for manipulation such as PCSK9, where identification of individuals with knockout mutations and benign lower LDL cholesterol has led to promising results in clinical trials and the development of evolocumab (Amgen) and alirocumab (Regeneron). In clinical practice and biomedical research, next-generation sequencing (NGS) and subsequent identification of genomic variants including single nucleotide variations, small insertions or deletions, and structural variants is an established method used to investigate the genetic causes and associations of disease. While whole genome and whole exome sequencing is a highly cost-effective and versatile method that assays gene sequence thus yielding both genetic and functional information. We, therefore, developed GeneTarget, a clinical genetic framework for the clinical trials analysis and interpretation of DNA sequencing data. The software is designed for use by clinicians and other users without an in-depth background in genetics.


## Methods
  
As depicted in Figure ??, VariantSurvival requires three types of input data: a neurological disease, treatment group meta data and an annotated multi-sample VCF file of the study group. As a first input, the user has to select from a list of neurological diseases. Following this selection, VariantSurvival curates a list of genes that are known to be associated with the selected disease. Additionally, VariantSurvival suggests SV calling methods that are beneficial to identify variants known to be present in the curated gene list. \
The second input is a set of all annotated structural variants from the study group. VariantSurvival requires a multi-sample VCF file (short “VCF file” from here) where each variant record was annotated with gene identifiers according to the Ensembl (ref) database. From the VCF file, a tally is created of all variant records that match the selected gene. Provided that the matching results in a non-empty set of SVs, the workflow continues with the group labeling.
The third input are the study group labels. Each patient within the tally is categorized with the study group label. The final labeled tally provides the input data for the survival analysis (ref). Here, VariantSurvival computes the cox regression (ref) to determine the difference in survival between the placebo and drug study group using the structural variant counts as covariates in the cox regression model. \
The green box of the pipeline in Figure 1 shows a generic predecessor workflow in order to create the second input for VariantSurvival. This predecessor workflow is not implemented in VariantSurvival. Supplementary section ?? provides recommendations and details about how to generate the required input formats for VariantSurvival.


###  Implementation
The VariantSurvival shiny app on GitHub https://github.com/collaborativebioinformatics/ VariantSurvival/ VariantSurvival_App.

The repository provides detailed instructions for tool usage and installation. 

The required packages to be installed before running the app are in the “requirement.txt” file.

The App require 1,5G of Ram and it is plateform independent.

The required inputs are : The variant VCF file and a metadata file. Examples of the files are in the “demo” folder.


###  Operation
The workflow of our tool is described as follows: As an initial step, the user will have the option to choose a disease from a list of neurological conditions ( Alzheimer's disease, Amyotrophic lateral sclerosis, Friedreich ataxia, Huntington's disease, Lewy body disease, Parkinson's disease, and Spinal muscular atrophy).
Once this is done, a list of genes known to be associated with the chosen disease is generated. 

The user needs to choose the target gene from the gene list. The structural variants count in the target gene region will be represented in a barplot. To verify the of the structural variants before starting the survival analysis.

Only the SVs that are in the target gene are considered significant, as a factor in the survival analysis. The placebo and treatment groups are identified using the metadata file.

The survival analysis result are represented in the second tab of the app. The first plot, compare the survival of the placebo and treatment group. The existence or not of the SVs is a factor, however, the count of the SVs is not considered.
The second plot illustrate the survival of the placebo and treatment group according to the SVs count in the target gene.
## Integrated Tools

Multiple R packages are combined to develop the shiny app including : shiny, shinydashboard, DT, vcfR, readr, readxl.

The survival analysis is performed and illustrated using the following R packages: survival, survminer, lubridate, gtsummary, ggsurvfit, dplyr, tidyverse, ggplot2.

The packages citation are in the "References.txt" file.
##  Flowchart

<img src="https://github.com/collaborativebioinformatics/GeneTarget/blob/main/img/GeneTargetWorkflow.svg" width="700">

# Shiny App

### Shiny app development

The apps interface was developed using multiple R packages including shiny (Winston Chang
et al., 2022), shinydashboard  (Winston Chang et al., 2021).  Multiple other R packages are intagrated including : DT (Yihui Xie et al, 20222), vcfR (Knaus BJ and  Grünwald NJ, 2017), readr (Wickham H et al 2022), readxl (Wickham H, and  Bryan J, 2022). .
The survival analysis
is performed and illustrated using the following R packages: survival (Therneau T, 2022; Terry M et al., 2000), survminer (Alboukadel Kassambara et al., 2021), lubridate (Grolemund G, and Wickham H, 2011), gtsummary(Sjoberg D et al., 2021), ggsurvfit (Daniel D. Sjoberg, 2022), dplyr (Wickham H et al., 2022), tidyverse (Wickham H et al., 2019), ggplot2 (Wickham H, 2016).

#### Shiny app interface

The workflow of our tool is described as follows: As an initial step, the user will have the option to choose a disease from a list of neurological conditions ( Alzheimer's disease, Amyotrophic lateral sclerosis, Friedreich ataxia, Huntington's disease, Lewy body disease, Parkinson's disease, and Spinal muscular atrophy).
Once this is done, a list of genes known to be associated with the chosen disease is generated. 

The user needs to choose the target gene from the gene list. The structural variants counted in the target gene region will be represented in a barplot to verify the presence of structural variants before starting the survival analysis.
From vcf to sv_event workflow    @Mostafa 
Only the SVs that are in the target gene are considered, as a factor in the survival analysis. The placebo and treatment groups are identified using the metadata file.

The survival analysis result are represented in the second tab of the app. The first plot, compare the survival of the placebo and treatment group. The existence or not of the SVs is a factor, however, the count of the SVs is not considered.
The second plot illustrate the survival of the placebo and treatment group according to the SVs count in the target gene.


## Input

### Shiny app interface

We present VariantSurvival, a lightweight shiny dashboard to visualize genotype-treatment response. The dashboard’s first tab provides all functionality to import VCF data, metadata and to choose a neurological disease as well as one or multiple target genes.


![image](https://user-images.githubusercontent.com/73958439/195391270-5c88d73c-1149-4fa5-af31-33e95d41a953.png)

#### Figure1. The VariantSurvival interface (First tab)

## Output


The second tab visualizes the survival analysis’ results:


![image](https://user-images.githubusercontent.com/73958439/195386048-2a978a86-9ba1-4834-b387-f320648de0fb.png)

#### Figure2. The VariantSurvival interface (Second tab)

# Team

* Ahmad Al Khleifat

* Thomas Krannich

* Hiba Ben Aribi

* Marina Herrera Sarrias

* Moustafa Shokrof
