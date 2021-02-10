#' @title P-value significance stars definer.
#'
#' @description Converts a p-value score in equivalent stars of significance.
#'
#' @param p.value A single numeric value indicating the p-value to evaluate.
#' @param one A numeric value to define the p-value threshold for the first level of significance (*). By default 0.05.
#' @param two A numeric value to define the p-value threshold for the second level of significance (**). By default 0.01.
#' @param three A numeric value to define the p-value threshold for the third level of significance (***). By default 0.001.
#' @param four A numeric value to define the p-value threshold for the fourth level of significance (****). By default 0.0001.
#'
#' @return It returns a string with the corresponding level of significance: \code{NS}, \code{*}, \code{**}, \code{***}, \code{****}.
#'
#' @examples
#' significance = pStars(0.002)
#'
#' require(dplyr)
#' data.frame =
#'    data.frame %>%
#'    mutate(p.stars = sapply(data.frame$p.value.column, pStars))
#'
#' @export pStars




pStars = function(p.value,
                  one = 0.05,
                  two = 0.01,
                  three = 0.001,
                  four = 0.0001) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  if (length(p.value) > 1) {return(warning("The 'p.value' parameter must a single numeric value."))}

  if (one > two & one > three & one > four &
      two > three & two > four &
      three > four) {
  if (p.value > one | is.na(p.value) | is.null(p.value)) {
    return("NS")} else {
      if (p.value <= one & p.value > two) {
        return("*")} else {
          if (p.value <= two & p.value > three) {
            return("**")} else {
              if (p.value <= three & p.value > four) {
                return("***")} else {
                  if (p.value <= four) {
                    return("****")}
                }
            }
        }
    }
  } else {return(warning("Limits must be in descending order, i.e. 0.05 > 0.01 > 0.001 > 0.0001"))}
}
