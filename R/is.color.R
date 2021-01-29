#' @title is.color
#' 
#' @description Function to define if each element of a string vector is an R-supported color string.
#' 
#' @param x A string vector.
#' 
#' @return A logical vector of the same length of \code{x}.
#' 
#' @export is.color



is.color = function(x) {
  
  # Check x class
  if (class(x) != "character") {return(warning("The input must be a string vector."))}
  
  return(sapply(x,
                function(X) {
                  tryCatch(is.matrix(col2rgb(X)), 
                           error = function(e) FALSE)},
                USE.NAMES = F))
}