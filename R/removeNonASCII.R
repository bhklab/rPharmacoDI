#' Remove all non-ascii characters from a data.frame
#'
#' @param df [`data.frame`] to purge of non-ASCII character
#'
#' @return [`data.frame`] with all character columns stripped of non-ASCII characters
#'
#' @import data.table
#' @export
removeNonASCII <- function(df) {
  DT <- data.table(df, keep.rownames='rn')
  charCols <- colnames(DT)[vapply(DT, is.character, logical(1))]
  DT[, (charCols) := lapply(.SD, iconv, from='ascii', to='ascii', sub=''), .SDcols=charCols]
  return(data.frame(DT[, -'rn'], row.names=DT$rn))
}