#' `get_disease_gene_list` 
#'
#' @param input_disease 
#' @return  genes_list 
get_disease_gene_list <- function(input_disease) {
  if(as.numeric(input_disease == 1)){
    genes_list  <-read_csv("disease_gene/ALS/genes_list.txt")
  }
  else if(as.numeric(input_disease == 2)){
    genes_list  <-read_csv("disease_gene/PD/genes_list.txt")
  }
  else if(as.numeric(input_disease == 3)) {
    genes_list  <-read_csv("disease_gene/AD/genes_list.txt")
  }
  else if(as.numeric(input_disease == 4)) {
    genes_list  <-read_csv("disease_gene/FD/genes_list.txt")
  }
  else {
    genes_list  <-read_csv("disease_gene/DLB/genes_list.txt")
  }
  return(genes_list)
}

#' `getGeneName` 
#'
#' @param info 
#' @return   
getGeneName <- function(info) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}

#' `CountSVsDf` 
#'
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
                       vcf){
  sample_disease_gene_df <- data.frame(
    matrix(0,
           ncol = ncol,
           nrow = nrow 
    )
  )
  colnames(sample_disease_gene_df) = disease_genes_names
  rownames(sample_disease_gene_df) = sample_names
  
  for (sv_idx in 1:length(genes_with_svs_in_sample)) {
    # if the sv is in a gene of interest
    sv_gene_name = genes_with_svs_in_sample[sv_idx]
    if(is.na(sv_gene_name))
      next;
    for(individual_idx in 1:length(samples)){
      indiv_id = samples[individual_idx]
      gt=vcf@gt[sv_idx, indiv_id]
      if(!is.na(gt))
      {
        sample_disease_gene_df[indiv_id, sv_gene_name] %+=% 1
      }
    }
  }
  sample_disease_gene_df <- rownames_to_column(sample_disease_gene_df,
                                               "patient_ID")
  return(sample_disease_gene_df)
}


#' `hist_df` 
#'
#' @param df 
#' @param input
#' @return  
hist_df <- function(df,input){
  svs_gene_input_df <- df[c("patient_ID",
                            input$target_gene,
                            input$phenotype)
                          ]
  colnames(svs_gene_input_df) <- c('patient_ID',
                                   'SV_count_per_gene',
                                   'Phenotype')
  
  svs_gene_input_df <- svs_gene_input_df %>%
    mutate(Phenotype = ifelse(Phenotype=="0",
                              "Placebo", 
                              "Treatment" )
           )
  return(svs_gene_input_df)
}

#' `RemoveNAs` 
#'
#' @param info 
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

#' implementation of += operator
#' https://stackoverflow.com/questions/5738831/
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

