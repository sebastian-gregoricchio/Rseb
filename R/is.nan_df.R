#' @title \code{is.nan()} applied to a data.frame
#'
#' @description Applies the function \code{is.nan()} to a full data.frame.
#'
#' @param data.frame Input data.frame.
#'
#' @return It returns a matrix/array containing logic values for each element of the input data.frame. When \code{TRUE} it means that the corresponding element is a \code{NaN}.
#'
#' @examples
#' is.nan.df(mtcars)
#'
#' @export is.nan_df


is.nan_df = function(data.frame) {
  do.call(cbind, lapply(data.frame, is.nan))}
