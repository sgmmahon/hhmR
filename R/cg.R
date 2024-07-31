#' @title cg
#'
#' @description Creates colour gradient between two hexcodes.
#'
#' @param colour1 The first hexcode colour.
#' @param colour2 The second hexcode colour.
#' @param n The length of the vector returned by the function.
#'
#' @return A vector of hexcodes of length n, containing a colour gradient between colour =1 and colour2.
#'
#' @examples
#' cg("white","black",20)
#'
#' @export
#' @importFrom grDevices colorRampPalette
cg = function(colour1, colour2, n = 15) {

  # Create a color palette function
  colour_func <- grDevices::colorRampPalette(c(colour1, colour2))

  # Generate the color gradient
  colour_gradient <- colour_func(n)

  # Return colour gradient
  return(colour_gradient)
}
