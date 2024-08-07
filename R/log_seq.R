utils::globalVariables(".")
#' @title log_seq
#'
#' @description Creates a vector of logarithmicly increasing values between 0 and a specified value `n`.
#' If `n` is specified as 1, the vector will be scaled to between 0 and 1.
#'
#' @param n The maximum value that the values in the sequence are scaled to.
#' @param ln How long the vector should be (defaults to 15).
#' @param round_values Option to round values to whole numbers (defaults to `TRUE`).
#' @param rmv_extremes Option to remove zero and the maximum value (i.e. `n`) from the beginning
#' and the end of the returned vector (defaults to `FALSE`). Note that this will mean the length
#' of the returned vector will be `n` - 2.
#'
#' @return A vector containing logarithmicly increasing values between 0 and a specified value `n`.
#'
#' @examples
#' # Create sequence of length 20, scaled between 0 and 500
#' log_seq(500,20)
#'
#' # Create sequence of length 15, scaled between 0 and 1
#' log_seq(1,12)
#'
#' @export
#' @importFrom utils globalVariables
log_seq = function(n,ln=15,round_values=T,rmv_extremes=F) {

  # Create a sequence from 1 to n.
  # If `n` is specified as 1, the vector will be scaled to between 0 and 1.
  if (n == 1) {
    seq = seq(1, 1000, length.out = ln)
    round_values = F
  } else {
    seq = seq(1, n, length.out = ln-1)
  }

  # Apply the logarithm to the sequence
  log_seq = log(seq)

  # Scale the sequence to the range [0, 1]
  min_val = min(log_seq)
  max_val = max(log_seq)
  one_seq = ((log_seq - min_val) / (max_val - min_val))

  # Reverse pattern of scale so breaks are focussed on lower rather than upper end of scale
  one_seq_rev = one_seq %>% {. - max(.)} %>% {. * -1} %>% rev()

  # Scale sequence to n
  if (round_values) {
    scaled_one_seq = one_seq_rev %>% {. * n} %>% round() %>% .[2:length(.)] %>% c(0,1,.)
  } else if (n == 1) {
    scaled_one_seq = one_seq_rev
  } else {
    scaled_one_seq = one_seq_rev %>% {. * n} %>% .[2:length(.)] %>% c(0,.Machine$double.xmin,.)
  }

  # Option to remove zero and the maximum value (i.e. `n`) from the beginning and the end of the vector
  if (rmv_extremes) {
    scaled_one_seq = scaled_one_seq %>% .[2:(length(.)-1)]
  }

  # Return scaled sequence
  return(scaled_one_seq)
}
