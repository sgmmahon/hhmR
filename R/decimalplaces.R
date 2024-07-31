#' @title decimalplaces
#'
#' @description Tests the number of non-zero decimal places within a number.
#'
#' @param x The number for the number of decimal places is to be measured.
#'
#' @return A single number, indicating the number of non-zero decimal places in `x`.
#'
#' @examples
#' decimalplaces(23.43234525)
#' decimalplaces(334.3410000000000000)
#' decimalplaces(2.000)
#'
#' @export
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
