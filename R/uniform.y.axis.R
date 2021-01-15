#' @title Plot Y-axis uniforming
#'  
#' @description Takes a list of ggplot2 plots, compares their Y-axis ranges and applies the highest/lowest limits to each plot in order to uniform all the plots. It can be used also to set the ticks step (to just change the breaks set all parameters as \code{FALSE}).
#' 
#' @param plot.list A single plot or a list of plots.
#' @param y.min Either a logical value to define whether uniform the lower limit or a numeric value defining the lower limit. By default \code{TRUE}.
#' @param y.max Either a logical value to define whether uniform the upper limit or a numeric value defining the upper limit. By default \code{TRUE}.
#' @param ticks.each Numeric value to define every how much should be placed a tick. By default \code{NULL}, ticks will be placed automatically.
#' 
#' @return Returns a plot list (or a single plot when only one input plot is provided) equivalent to the input list provided by the user in which the Y-axis of all the plots will be uniformed.
#' 
#' @export uniform.y.axis
 
uniform.y.axis = function(
  plot.list,
  y.min = TRUE,
  y.max = TRUE,
  ticks.each = NULL) 
  
{ # BEGIN function -----------------------------------------------------------------------
  
  # Check list length
  if (!("list" %in% class(plot.list))) {
    if ("ggplot" %in% class(plot.list) | "gg" %in% class(plot.list)) {
      plot.list = list(plot.list)
    }  else {
      return(warning("The 'plot.list' parameter must be a ggplot object or a list of ggplot objects."))}
  }
  
  
  
  # Load required libraries
  require(ggplot2)
  
  
  
  # Get current Y ranges
  y.min_list = c()
  y.max_list = c()
  
  for (i in 1:length(plot.list)) {
    y.min_list = c(y.min_list, ggplot_build(plot.list[[i]])$layout$panel_params[[1]]$y$limits[1])
    y.max_list = c(y.max_list, ggplot_build(plot.list[[i]])$layout$panel_params[[1]]$y$limits[2])
  }
  
  
  
  # Re-define the limits
  if (length(y.min) != 1 | length(y.max) != 1) {return(warning("The 'y.min' and 'y.max' parameters must be either a unique logical or numeric value."))}
  
  if (is.na(y.min) | is.numeric(y.min) | is.null(y.min)) {
    new_y.min = y.min} else if (y.min == TRUE) {
      new_y.min = min(y.min_list)} else {
        new_y.min = NA}
  
  if (is.na(y.max) | is.numeric(y.max) | is.null(y.max)) {
    new_y.max = y.max} else if (y.max == TRUE) {
      new_y.max = max(y.max_list)} else {
        new_y.max = NA}

  
  
  # Re-define breaks
  if (!is.null(ticks.each)) {
    breaks_list = list()
    
    for (i in 1:length(plot.list)) {
      min = ifelse(test = is.na(new_y.min), yes = y.min_list[i], no = new_y.min)
      max = ifelse(test = is.na(new_y.max), yes = y.max_list[i], no = new_y.max)
      breaks_list[[i]] = seq(min, max, ticks.each)
    }
  }
  

  
  # Re-define the Y-axis
  for (i in 1:length(plot.list)) {
    if (is.null(ticks.each)) {
      plot.list[[i]] = plot.list[[i]] + scale_y_continuous(limits = c(new_y.min, new_y.max))
    } else {
      plot.list[[i]] = plot.list[[i]] + scale_y_continuous(breaks = breaks_list[[i]], limits = c(new_y.min, new_y.max))
    }
  }
  
  
  # Return the modified plots
  if (length(plot.list) == 1) {return(plot.list[[1]])} else {return(plot.list)}

  
} # END function
