#' map_col_names
#'
#' @param input input
#' @param cov_list list
#' @noRd
#' @return cov_list
map_col_names <- function(input, cov_list){
  if(input$event %in% cov_list){
    cov_list <- str_replace(cov_list,
                            input$event,
                            "event")}
  if(input$group %in% cov_list){
    cov_list <- str_replace(cov_list,
                            input$group,
                            "trial_group_bin")
  }
  return(cov_list)
}
