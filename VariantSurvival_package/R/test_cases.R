############# Helper functions ###################
install_load_requirements<- function() {
  if (!require("shiny")) install.packages("shiny")
  if (!require("shinyjs")) install.packages("shinyjs")
  if (!require("shinyWidgets")) install.packages("shinyWidgets")
  if (!require("shinythemes")) install.packages("shinythemes")
  if (!require("shinycssloaders")) install.packages("shinycssloaders")
  if (!require("shinydashboard")) install.packages("shinydashboard")
  if (!require("DT")) install.packages("DT")
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("vcfR")) install.packages("vcfR")
  if (!require("readr")) install.packages("readr")
  if (!require("ggsurvfit")) install.packages("ggsurvfit")
  if (!require("survival")) install.packages("survival")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("survminer")) install.packages("survminer")
  if (!require("lubridate")) install.packages("lubridate")
  if (!require("gtsummary")) install.packages("gtsummary")
  if (!require("tools")) install.packages("tools")
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(shinythemes)
  library(shinyWidgets)
  library(shinycssloaders)
  library(DT)
  library(dplyr)
  library(tidyverse)
  library(vcfR)
  library(readr)
  library(readxl)
  library(ggsurvfit)
  library(survival)
  library(ggplot2)
  library(survminer)
  library(lubridate)
  library(gtsummary)
}


#' `get_disease_gene_list` parses an existing .txt
#' file containing all known target genes associated
#' with the selected disease.
#' @param input_disease: input disease selection
#' @return  genes_list
get_disease_gene_list <- function(disease_gene, input_disease) {
  genes_list <-( na.omit(disease_gene[,input_disease])
                 %>% rename(genes = input_disease))
  return(genes_list$genes)
}


#' `checkInput`
#' @param info
#' @return
checkInput <- function(input) {
  if (input$ids != "N/A"
      & input$time != "N/A"
      & input$group != "N/A"
      & input$event != "N/A"
      & input$target_gene != "N/A"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#' `getGeneName`
#' @param info
#' @return
getGeneName <- function(info, geneIDS) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}


get_cov_list <- function(metadata, input) {
  inputs_list <- c(input$event, input$time, input$ids)
  cov_list <- colnames(metadata)
  cov_list <- cov_list[!(cov_list %in% inputs_list)]
  return(cov_list)
}


#' `CountSVsDf`
#' @param ncol
#' @param nrow
#' @param disease_genes_names
#' @param  sample_names
#' @param genes_with_svs_in_sample
#' @param vcf
#' @return
CountSVsDf <- function(ncol, nrow,
                       disease_genes_names,
                       sample_names,
                       genes_with_svs_in_sample,
                       vcf,
                       ids_col){
  sample_disease_gene_df <- data.frame(matrix(0, ncol = ncol, nrow = nrow))
  colnames(sample_disease_gene_df) <- disease_genes_names
  rownames(sample_disease_gene_df) <- sample_names
  
  for (sv_idx in 1:length(genes_with_svs_in_sample)) {
    # if the sv is in a gene of interest
    sv_gene_name <- genes_with_svs_in_sample[sv_idx]
    if(is.na(sv_gene_name))
      next
    for(individual_idx in 1:length(sample_names)){
      indiv_id <- sample_names[individual_idx]
      gt <- vcf@gt[sv_idx, indiv_id]
      if(!is.na(gt))
      {
        sample_disease_gene_df[indiv_id, sv_gene_name] %+=% 1
      }
    }
  }
  sample_disease_gene_df <- (rownames_to_column(sample_disease_gene_df,ids_col)
                             %>%  rename(ids = ids_col))
  return(sample_disease_gene_df)
}


#' `RemoveNAs` removes all rows with Na
#' entries in the time factor column.
#' @param df: metadata data frame
#' @param time_col: time factor input
#' @return
RemoveNAs <- function(df, time_col) {
  entries <- df[[time_col]]
  nas_entries <- is.na(entries)
  NAs_entries <- entries=="NA"
  if (sum(nas_entries) > 0){
    return(df[!nas_entries,])
  } else if (sum(NAs_entries) > 0){
    return(df[!NAs_entries,])
  } else {
    return(df)
  }
}

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

#' implementation of += operator
#' https://stackoverflow.com/questions/5738831/
`%+=%` <- function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

############# Test cases ###################

# load requirements
install_load_requirements()

# input file requirements: 
# ensembleTogenes.csv
# merged.filtered.vcf
# metadata.xlsx
######################################################
x <- "merged.filtered.vcf"
y <- "metadata.xlsx"
z <- 'ensembleTogenes.csv'
metadata <- "metadata.xlsx"
input.patient_id <- "patient_ID"
input.time <- "Time_to_death_or_last_followup_years"
input.surv_status_bin <-  "survival.status_bin"
input.phenotype <- "Phenotype"
input.gene_id <- "SETX"
#######################################################################################

vcf <- vcfR::read.vcfR(x, verbose = FALSE)
metadata <- readxl::read_excel(metadata)
metadata <- na.omit(metadata)
# VCF genotype information:
sample_names = colnames(vcf@gt)[-1] # exclude the first element "FORMAT"
gene_ids_table <- read.csv(file = z)
rownames(gene_ids_table) <- gene_ids_table$ensembleID
gene_names = gene_ids_table$GeneName
gene_ids = gene_ids_table$ensembleID
disease_genes_names = gene_ids_table$GeneName
genes_with_svs_in_sample <- apply(vcf@fix, 1, getGeneName, gene_ids_table)
count_df <- CountSVsDf(length(disease_genes_names),
                       length(sample_names),
                       disease_genes_names,
                       sample_names,
                       genes_with_svs_in_sample,
                       vcf,
                       input.patient_id)
new_md <- (RemoveNAs(metadata, input.time)
           %>% rename(ids = patient_ID,
                      trial_group_bin = input.phenotype,
                      event = input.surv_status_bin,
                      time = input.time)
           )
new_df <- (subset(count_df, ids %in% new_md$ids)
           %>% rename(SV_count_per_gene = input.gene_id)
           )
no_na_df <- merge(new_df, new_md, by = "ids")
no_na_df <- transform(no_na_df,
                      time = as.numeric(time),
                      event = as.numeric(event),
                      SV_count_per_gene = as.numeric(SV_count_per_gene)
                      )
no_na_df <- (no_na_df
             %>% mutate(SV_bin = ifelse(new_df$SV_count_per_gene>0, 1, 0))
             %>% mutate(SV = ifelse(SV_bin=="0", "with", "without"))
             %>% mutate(trial_group = ifelse(trial_group_bin=="0",
                                             "Placebo",
                                             "Treatment" )
                        )
             )
count_tab <- as.data.frame(apply(new_df[,-1], 2, function(c)sum(c!=0)))
colnames(count_tab) <- "Count of patients with SVs"
count_tab <- rownames_to_column(count_tab, "Gene ID")
count_tab <- count_tab %>% arrange(desc(count_tab['Count of patients with SVs']))
#####################################################################################
## KM
#####################################################################################

#-------------------------Null model--------------------------------------------------

null_model <- survfit(Surv(time, event)~1, data = no_na_df, type = "kaplan-meier")

null_model_plot <- ggsurvplot_combine(list(null_model),
                                      data=no_na_df,
                                      risk.table = TRUE,
                                      conf.int = FALSE,
                                      conf.int.style = "step",
                                      risk.table.y.text = FALSE,
                                      xlab = "Years",
                                      ggtheme = theme(
                                        text = element_text(size = 20),
                                        panel.background = element_rect(fill = "white",
                                                                        colour = "white"),
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major.y = element_line(colour='white'),
                                        panel.grid.minor.y = element_line(colour='white')
                                      )
)
# life table
lt_null_model <- round_df(as.data.frame(surv_summary(null_model)), 3)


#-------------------------Multiple model--------------------------------------------------#

sc_without <- survfit2(Surv(time, event)~trial_group_bin, data = no_na_df, 
                       subset=(SV_bin==0), type = "kaplan-meier")
sv_with <- survfit2(Surv(time, event)~trial_group_bin,data = no_na_df,
                    subset=(SV_bin==1), type = "kaplan-meier")
surv_fit_list <- list("with SV" = sv_with, "without SV" = sc_without)
groups_df <- (unique(no_na_df[c("SV_bin", "trial_group_bin")])
              %>% arrange(SV_bin, trial_group_bin))
lables <- c()
cols <- c()
for (row in 1:nrow(groups_df)){
  sv_bin <- groups_df[row, "SV_bin"]
  trial_g <- groups_df[row, "trial_group_bin"]
  if (sv_bin == 0){
    if (trial_g == 0){
      lables <- c(lables, "without sv - placebo")
      cols <- c(cols, "Violetred2")
    }
    else if (trial_g == 1){
      lables <- c(lables, "without sv - treatment")
      cols <- c(cols, "turquoise3")
    }
  }
  if (sv_bin == 1){
    if (trial_g == 0){
      lables <- c(lables, "with sv - placebo")
      cols <- c(cols, "Violetred4")
    }
    else if (trial_g == 1){
      lables <- c(lables, "with sv - treatment")
      cols <- c(cols, "steelblue")
    }
  }
}

mm_model_plot <- ggsurvplot_combine(surv_fit_list,
                                    data=no_na_df,
                                    risk.table = TRUE,
                                    conf.int = FALSE,
                                    conf.int.style = "step",
                                    risk.table.y.text = FALSE,
                                    xlab = "Years",
                                    ggtheme = theme(
                                      text = element_text(size = 20),
                                      panel.background = element_rect(fill = "white",
                                                                      colour = "white"),
                                      axis.line = element_line(colour = "black"),
                                      panel.grid.major.y = element_line(colour='white'),
                                      panel.grid.minor.y = element_line(colour='white'),
                                      # legend.position = c(0.2, 0.5)
                                    ),
                                    legend = "right",
                                    legend.title = "Group",
                                    legend.labs =lables,
                                    palette = cols)

# life tables:
tl_sv_without <- round_df(as.data.frame(surv_summary(sc_without, data = no_na_df)), 3)
tl_sv_with <- data.frame(lapply(surv_summary(sv_with, data = no_na_df),
                                function(x) if(is.numeric(x)) round(x, 3) else x))

tl_sv_with <- tl_sv_with[, !names(tl_sv_with) %in% c("strata")]

