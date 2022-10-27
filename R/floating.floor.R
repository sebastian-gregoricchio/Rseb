#' @title Flooring to floating values
#'
#' @description Computes the floor of the given value but with any number of digits (to the closest floating number of given digits).
#'
#' @param num A single number or a numeric vector.
#' @param digits A single integer indicating the maximum number of digits required.
#'
#' @return A floored number or numeric vector.
#'
#' @export floating.floor

floating.floor = function(num,
                          digits = 1) {

  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)

  # Check that digits is an integer
  if (digits != round(digits)) {return(warning("The 'digits' parameter must be an integer."))}

  # Avoiding negative digits and calculating the number length
  digits = nchar(as.character(trunc(abs(num)))) + abs(digits)

  # Considering the abs(numbers) < 1
  if (abs(num) < 1) {
    num = num * 10
    digits = digits - 1
  }

  # Handle for digits 0 digits, equivalent of normal floor
  power = floor(log10(abs(num))) + 1 - digits
  rounded.num = floor(num / 10^power) * 10^power

  # Handle the num = 0 case
  rounded.num[num == 0] = 0

  if (abs(num/10) < 1) {
    return(rounded.num/10)}
  else {
    return(rounded.num)}

} # End function
