utils::globalVariables(".")
#' @title exp_seq
#'
#' @description Creates a vector of exponentially increasing values between 0 and a specified value `n`.
#' If `n` is specified as 1, the vector will be scaled to between 0 and 1.
#'
#' @param n The maximum value that the values in the sequence are scaled to.
#' @param ln How long the vector should be (defaults to 15).
#' @param exponent The exponential power with which to multiply the sequence by (defaults to 2).
#' @param round_values Option to round values to whole numbers (defaults to `TRUE`). If `n` equals 1,
#' round_values will automatically be set to FALSE.
#' @param rmv_extremes Option to remove zero and the maximum value (i.e. `n`) from the beginning
#' and the end of the returned vector (defaults to `FALSE`). Note that this will mean the length
#' of the returned vector will be `n` - 2.
#'
#' @return A vector containing exponentially increasing values between 0 and a specified value `n`.
#'
#' @examples
#' # Create sequence of length 8, scaled between 0 and 10000
#' exp_seq(10000,8)
#' # Set rmv_extremes = FALSE to get full sequence
#' exp_seq(10000,8,rmv_extremes = FALSE)
#' # The exponent defaults to 2. Setting it to between 1 and 2 causes it to converge on
#' # a linear sequence. When exponent is set to 1 the sequence increases linearly
#' exp_seq(10000,8,exponent=1)
#' # Setting it to greater than 2 will cause it the values in the sequence to shift towards zero
#' exp_seq(10000,8,exponent=4)
#'
#' # Create sequence of length 12, scaled between 0 and 1
#' exp_seq(1,12)
#' exp_seq(1,12,rmv_extremes = FALSE)
#' exp_seq(1,12,exponent=1)
#' exp_seq(1,12,exponent=4)
#'
#' @export
#' @importFrom utils globalVariables
exp_seq = function(n,ln=15,exponent=2,round_values=TRUE,rmv_extremes=TRUE) {

  # Check variables have been entered correctly
  if (!is.numeric(exponent) || exponent < 1) { stop("`exponent` must be a numer greater than or equal to 1.") }

  # Create a sequence from 1 to n.
  # If `n` is specified as 1, the vector will be scaled to between 0 and 1.
  if (n == 1) {
    seq = seq(0, 1000, length.out = ln)
    round_values = FALSE
  } else {
    seq = seq(0, n, length.out = ln)
  }

  # Apply the logarithm to the sequence
  exp_seq = seq %>% {.^exponent}

  # Scale the sequence to the range [0, 1]
  min_val = min(exp_seq)
  max_val = max(exp_seq)
  one_seq = ((exp_seq - min_val) / (max_val - min_val))

  # Scale sequence to n
  if (round_values) {
    scaled_one_seq = one_seq %>% {. * n} %>% round()
  } else if (n == 1) {
    scaled_one_seq = one_seq
  } else {
    scaled_one_seq = one_seq %>% {. * n}
  }

  # Option to remove zero and the maximum value (i.e. `n`) from the beginning and the end of the vector
  if (rmv_extremes) {
    scaled_one_seq = scaled_one_seq %>% .[2:(length(.)-1)]
  }

  # Return scaled sequence
  return(scaled_one_seq)
}
