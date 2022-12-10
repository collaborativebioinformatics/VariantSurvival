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
