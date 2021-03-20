#' @title Plot density signal of NGS data.
#'
#' @description Plots the density profile of NGS data (e.g. ChIP-seq, ATAC-seq, MeDIP-seq, etc.). Used by the function \code{\link{plot.density.profile}}.
#'
#' @param samples A character vector containing the samples list.
#' @param scores A numeric vector containing the scores for the Y-axis.
#' @param positions A numeric vector containing the position for the X-axis.
#' @param variance_scores A numeric vector containing the variance/error value at each position.
#' @param xlab A string containing the label for the X-axis. By default "Distance from regions center [bp]".
#' @param ylab A string containing the label for the Y-axis. By default "Average density signal".
#' @param line_type Vector to define each line type. Both numeric and string codes are accepted. if only one element is given this will be applied to all the lines. By default "solid". \cr Example 1: \code{c("solid", "dashed")}. \cr Example 2: \code{c(1, 2)}
#' @param x_lim List of numeric vectors with two elements each to define the range of the X-axis. To set only one side use NA for the side to leave automatic. If only one range is given this one will be applied to all the plots. By default \code{NULL}, the range will be defined automatically. \cr Example \code{list(c(0, 20), c(NA, 30), c(0, NA), c(NA, NA))}.,
#' @param y_lim List of numeric vectors with two elements each to define the range of the Y-axis. To set only one side use NA for the side to leave automatic. If only one range is given this one will be applied to all the plots. By default \code{NULL}, the range will be defined automatically. \cr Example \code{list(c(0, 20), c(NA, 30), c(0, NA), c(NA, NA))}.,
#' @param x_intercept A vector indicating the X intercepts for the vertical lines. By default 0.
#' @param colors Vector to define the line and error area colors. If only one value is provided or the number of values is lower than the required ones only the first value will be used. All standard R.colors values are accepted. By default \code{c("blue", "red", "purple", "orange", "green")}.
#' @param title A string containing the label for the X-axis. By default "Density profile".
#' @param text_size Numeric value to define the size of the text for the labels of all the plots. By default 12.
#' @param variance Logic value to define whether to plot the error/variance around the signal. By default \code{TRUE}.
#' @param print_plot Logic value to define whether to print the plot once generated or not. By default \code{FALSE}.
#' @param line_width Numeric value to define the line width for all the plots. By default 1.,
#' @param variance_opacity Numeric value to define the alpha/transparency of the error/variance. By default 0.25. Parameter considered only when \code{variance = TRUE)}.
#'
#' @return Returns a plot in ggplot2 format.
#'
#' @export density_plot
#'
# @import ggplot2


density_plot = function(
  samples,
  scores,
  positions,
  variance_scores,
  xlab = "Distance from regions center [bp]",
  ylab = "Average density signal",
  line_type = "solid",
  y_lim = NULL,
  x_lim = NULL,
  x_intercept = 0,
  colors = c("blue", "red", "purple", "orange", "green"),
  title = "Density profile",
  text_size = 12,
  variance = T,
  print_plot = F,
  line_width = 1,
  variance_opacity = 0.25) { # BEGIN

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#


  # Create a matrix
  matrix =
    data.frame(sample = samples,
               position = positions,
               score = scores,
               variance = variance_scores)

  # Set to 0 the variance if variance == F
  if (variance == F) {matrix$variance = 0}

  # scaling size of parameters
  n_samples = length(unique(matrix$sample))

  # scaling number of colors/line.type to be equal to number samples
  scaling = function(x, n_samples) {
    if (length(x) < n_samples &
        length(x) == 1) {
      v = c(rep(x, n_samples))
    } else if (length(x) > 1 &
               length(x) < n_samples) {
      v = c(rep(x[1], n_samples))
      message("Number of colors and/or line.type lower than number of samples --> The first value is applied to all samples/groups.")
    } else (v = as.vector(x))
    return(v)
  }

  line_type = scaling(line_type, n_samples)
  colors = scaling(colors, n_samples)

  require(ggplot2)

  # building the plot
  plot =
    ggplot(data = matrix,
           aes(x = position,
               y = score,
               ymin = score-variance,
               ymax = score+variance,
               group = sample,
               color = sample,
               fill = sample,
               linetype = sample)) +
    geom_line(size = line_width) +
    scale_color_manual(values = as.vector(colors[1:n_samples])) +
    scale_linetype_manual(values = as.vector(line_type[1:n_samples])) +
    scale_fill_manual(values = as.vector(colors[1:n_samples])) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    geom_vline(xintercept = x_intercept,
               linetype = "dashed",
               col = "gray70") +
    theme_classic() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(color = "#000000"),
          axis.title = element_text(color = "#000000"),
          axis.line = element_line(color = "#000000"),
          axis.ticks = element_line(color = "#000000"),
          title = element_text(color = "#000000"),
          legend.title = element_text(color = "#000000"),
          legend.text = element_text(color = "#000000"))

  # Add variance if required
  if (variance == T) {plot = plot + geom_ribbon(alpha = variance_opacity, color = NA)}


  if (!(is.null(y_lim))) plot = plot + ylim(y_lim)
  if (!(is.null(x_lim))) plot = plot + xlim(x_lim)

  if (print_plot == T) {
    (print(plot))}

  return(plot)

} # END
