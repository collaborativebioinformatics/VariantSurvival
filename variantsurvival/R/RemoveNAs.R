#' `RemoveNAs` removes all rows with Na
#' entries in the time factor column.
#' @param df: metadata data frame
#' @param time_col: time factor input
#' @return df
#' @noRd
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
