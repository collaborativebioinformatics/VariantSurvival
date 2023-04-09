#' get_cov_list
#'
#' @param metadata file
#' @param input file
#'
#' @return cov_list

get_cov_list <- function(metadata, input) {
  inputs_list <- c(input$event, input$time, input$ids)
  cov_list <- colnames(metadata)
  cov_list <- cov_list[!(cov_list %in% inputs_list)]
  return(cov_list)
}
