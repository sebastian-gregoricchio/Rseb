#' @title Data frame conversion to a list of columns.
#'
#' @description Converts each column of a data.frame in a element of a list with the corresponding name of the original column. Useful for further use in functions such as purrr::pmap().
#'
#' @param x A data.frame to be converted
#'
#' @return A list of vectors in which each element is a column of input the data.frame.
#'
#' @examples
#' data.frame.to.list(mtcars)
#'
#' @export data.frame.to.list


data.frame.to.list = function(x) {
  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  if (class(x) != "data.frame") {
    return(warning("The input must be a data.frame object"))}

  data_frame_list = lapply(X = c(1:ncol(x)), FUN = function(i)(return((x[,i]))))
  names(data_frame_list) = names(x)
  return(data_frame_list)
}
