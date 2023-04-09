#' `CountSVsDf`
#' @param ncol name column
#' @param nrow name row
#' @param disease_genes_names disease name
#' @param sample_names sample name
#' @param genes_with_svs_in_sample gene
#' @param vcf file
#' @return sample_disease_gene_df
#' @noRd
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
