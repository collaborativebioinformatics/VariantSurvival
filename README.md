# GeneTarget : A tool to identify genotype-treatment response

# Team
Ahmad/
Thomas/
Marina/
Hiba/
Moustafa/


##  Abstract
For a number of neurological diseases, such as Alzheimer’s disease, Parkinson’s disease and many others, there are genes known to be implicated in the mechanism of the disease. A common question in this matter is whether a novel structural variant in any such gene may be related to drug response, and how this information could contribute in the analysis of clinical trials. To this end, we introduce GeneTarget, a tool that will make it possible to identify gene treatment response using genotype data


##  Implementation
GeneTarget is available on GitHub https://github.com/collaborativebioinformatics/directed_sv_annotation. The repository provides detailed instructions for tool usage and installation. A bash script for an automated installation of the required dependencies.


##  Operation
The workflow of our tool is described as follows:
As an initial step, the user will have the option to choose a disease from a list of neurological conditions. Once this is done, a list of genes known to be associated with the chosen disease is generated. In parallel,  NGs are also required as inputs in order to generate a BAM alignment file, which is then used to generate a VCF file using a specific SV detection method. Subsequently, a filtering step is included, so only the SVs that are in the genes list are included. Provided that this intersection shows results, this step triggers the downstream pipeline.
The downstream pipeline starts with the labelling of clinical trial data, where a sample of individuals is selected and classified according to a treatment group (placebo/treatment). Next, tabular data regarding structural variant events are collected for all individuals. The latter are then used as input for the survival analysis step, which is performed on each of the treatment groups. 

##  Flowchart

![]![Screenshot 2022-10-10 at 18 42 54](https://user-images.githubusercontent.com/41301333/194926460-94f62ffd-71e3-48e5-a764-b28f57c69fac.png)

# Shiny App

![]![image (2)](https://user-images.githubusercontent.com/41301333/194929191-437bfa99-9a31-4b63-be02-4a4a874ceb29.png)


![image (5)](https://user-images.githubusercontent.com/41301333/194935518-d790a587-f328-4785-8156-6f1d8d6524f1.png)

