#' @title Gradient colors generation and assignment
#'
#' @description Give a vector of colors generates a finite number of shadows that will be assigned to a numeric vector depending on the value of each element.
#'
#' @param values A numeric vector containing the values to which a color must be assigned (NAs and NaN will be converted to 0).
#' @param colors A string vector with the colors, in the wished order, that have to be used to generated the shadows. By default \code{c("blue", "white", "red")}.
#' @param bins An atomic integer value to define the total number of bins/steps in which the gradient should be dived.
#'
#' @return A vector containing the assigned colors corresponding to each element of \code{values}.
#'
#' @export color.gradient


color.gradient = function(values,
                          colors = c("blue", "white", "red"),
                          bins = 100) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # Check values class
  if (class(values) %in% c("numeric", "integer")) {
    if (is.character(values)) {return(warning("The 'values' input must be a numeric vector."))}
  } else {return(warning("The 'values' input must be a numeric vector."))}


  # Check colors
  if (class(colors) != "character") {return(warning("The 'colors' parameter must be a vector of strings containing colors values in any R-supported format."))}
  if (F %in% Rseb::is.color(colors)) {return(warning("The 'colors' parameter must be a vector of strings containing colors values in any R-supported format."))}

  # Check steps
  if (!is.atomic(bins) | !is.numeric(bins)) {return(warning("The 'bins' value must be a single integer."))}


  # Convert NA and NaN values to 0
  if (TRUE %in% (is.na(values) | is.nan(values))) {message("The values provided contain NAs and/or NaNs that will be converted to 0.")}

  values = sapply(values, function(x){return(ifelse(test = is.na(x) | is.nan(x),
                                                    yes = 0,
                                                    no = x))})


  return(colorRampPalette(colors)(bins)[findInterval(values, seq(min(values), max(values), length.out = bins))])

} # END function
