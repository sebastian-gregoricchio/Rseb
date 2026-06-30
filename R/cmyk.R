#' @title CMYK color converter
#'
#' @description Converts CMYK color values to hexadecimal color values
#'
#' @param C Value in the 0-100 range for the Cyan component.
#' @param M Value in the 0-100 range for the Magenta component.
#' @param Y Value in the 0-100 range for the Yellow component.
#' @param K Value in the 0-100 range for the Key component.
#'
#' @return The result is a string for the color in hexadecimal scale, eg. "#FFFFFF".
#'
#' @examples
#' color = cmyk(0, 0, 0, 0)
#'
#' @export cmyk


cmyk = function(C, M, Y, K) {

  C = C / 100.0
  M = M / 100.0
  Y = Y / 100.0
  K = K / 100.0

  n.c = (C * (1-K) + K)
  n.m = (M * (1-K) + K)
  n.y = (Y * (1-K) + K)

  r.col = ceiling(255 * (1-n.c))
  g.col = ceiling(255 * (1-n.m))
  b.col = ceiling(255 * (1-n.y))

  return(rgb(r.col/255, g.col/255, b.col/255))
}
