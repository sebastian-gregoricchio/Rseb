#' @title Mode calculation
#'
#' @description Calculate the mode value of a vector of numeric values.
#'
#' @param v A vector of numeric numbers
#'
#' @return A single number corresponding to the mode of the list of numbers give as input
#'
#' @examples
#' mode = calculate.mode(v = c(6, 8, 4, 845, 8, 5, 55, 84, 8, 84, 45, 5))
#'
#' @export calculate.mode

calculate.mode = function(v) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  uniq.v = unique(v)
  uniq.v[which.max(tabulate(match(v, uniq.v)))]
}
