#' round_df
#' round all numeric variables
#' @param x data frame
#' @param digits number of digits to round
#' @noRd
#' @return x

round_df <- function(x, digits) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
