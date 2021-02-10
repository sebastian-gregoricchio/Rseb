#' @title Plot X-axis uniforming
#'
#' @description Takes a list of ggplot2 plots, compares their X-axis ranges and applies the highest/lowest limits to each plot in order to uniform all the plots. It can be used also to set the ticks step (to just change the breaks set all parameters as \code{FALSE}).
#'
#' @param plot.list A single plot or a list of plots.
#' @param x.min Either a logical value to define whether uniform the lower limit or a numeric value defining the lower limit. By default \code{TRUE}.
#' @param x.max Either a logical value to define whether uniform the upper limit or a numeric value defining the upper limit. By default \code{TRUE}.
#' @param ticks.each Numeric value to define every how much should be placed a tick. By default \code{NULL}, ticks will be placed automatically.
#' @param digits A single integer indicating the maximum number of digits required for the rounding of the axis values.
#'
#' @return Returns a plot list (or a single plot when only one input plot is provided) equivalent to the input list provided by the user in which the X-axis of all the plots will be uniformed.
#'
#' @export uniform.x.axis

uniform.x.axis = function(
  plot.list,
  x.min = TRUE,
  x.max = TRUE,
  ticks.each = NULL)

{ # BEGIN function -----------------------------------------------------------------------

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # Check list length
  if (!("list" %in% class(plot.list))) {
    if ("ggplot" %in% class(plot.list) | "gg" %in% class(plot.list)) {
      plot.list = list(plot.list)
    }  else {
      return(warning("The 'plot.list' parameter must be a ggplot object or a list of ggplot objects."))}
  }



  # Load required libraries
  require(ggplot2)



  # Get current X ranges
  x.min_list = c()
  x.max_list = c()

  for (i in 1:length(plot.list)) {
    x.min_list = c(x.min_list, ggplot_build(plot.list[[i]])$layout$panel_params[[1]]$x$limits[1])
    x.max_list = c(x.max_list, ggplot_build(plot.list[[i]])$layout$panel_params[[1]]$x$limits[2])
  }



  # Re-define the limits
  if (length(x.min) != 1 | length(x.max) != 1) {return(warning("The 'x.min' and 'x.max' parameters must be either a unique logical or numeric value."))}

  if (is.na(x.min) | is.numeric(x.min) | is.null(x.min)) {
    new_x.min = x.min} else if (x.min == TRUE) {
      new_x.min = min(x.min_list)} else {
        new_x.min = NA}

  if (is.na(x.max) | is.numeric(x.max) | is.null(x.max)) {
    new_x.max = x.max} else if (x.max == TRUE) {
      new_x.max = max(x.max_list)} else {
        new_x.max = NA}

  # Rund the new range
  if (!is.na(new_x.min)) {new_x.min = Rseb::floating.floor(new_x.min, digits = digits)}
  if (!is.na(new_x.max)) {new_x.max = Rseb::floating.ceiling(new_x.max, digits = digits)}


  # Re-define breaks
  if (!is.null(ticks.each)) {
    breaks_list = list()

    for (i in 1:length(plot.list)) {
      min = ifelse(test = is.na(new_x.min), yes = x.min_list[i], no = new_x.min)
      max = ifelse(test = is.na(new_x.max), yes = x.max_list[i], no = new_x.max)
      breaks_list[[i]] = seq(min, max, ticks.each)
    }
  }



  # Re-define the X-axis
  for (i in 1:length(plot.list)) {
    if (is.null(ticks.each)) {
      plot.list[[i]] = plot.list[[i]] + scale_x_continuous(limits = c(new_x.min, new_x.max))
    } else {
      plot.list[[i]] = plot.list[[i]] + scale_x_continuous(breaks = breaks_list[[i]], limits = c(new_x.min, new_x.max))
    }
  }


  # Return the modified plots
  if (length(plot.list) == 1) {return(plot.list[[1]])} else {return(plot.list)}

} # END function
