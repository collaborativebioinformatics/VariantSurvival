# Define server logic required to draw a histogram

#' Sum of vector elements
#'
#' `sum` returns the sum of all the values present in its arguments.
#' @param input 
#' @param output 
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


getGeneName <- function(info) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}

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


getGeneName <- function(info) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}


`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))