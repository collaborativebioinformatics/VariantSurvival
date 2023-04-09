#' genesCountTable
#'
#' @param vcf file
#' @param metadata file
#' @param input input
#' @param gene_ids_table table
#' @noRd
#' @return count_table
genesCountTable <- function(vcf, metadata, input, gene_ids_table){
  genes_with_svs_in_sample <- apply(vcf@fix, 1, getGeneName, gene_ids_table)
  sample_names <- colnames(vcf@gt)[-1]
  disease_genes_names <- gene_ids_table$GeneName
  count_df <- CountSVsDf(length(disease_genes_names),
                         length(sample_names),
                         disease_genes_names,
                         sample_names,
                         genes_with_svs_in_sample,
                         vcf,
                         input$ids)
  new_md <- RemoveNAs(metadata, input$time) %>% rename(ids = input$ids)
  # what happen if the input$target_gene is not in the vcf file?
  new_df <- subset(count_df, ids %in% new_md$ids)
  ##
  count_table <- as.data.frame(apply(new_df[,-1], 2, function(c)sum(c!=0)))
  colnames(count_table) <- "Count of patients with SVs"
  count_table <- rownames_to_column(count_table, "Gene_ID")
  count_table <- count_table %>% arrange(desc(count_table["Count of patients with SVs"]))
  return(count_table)
}
