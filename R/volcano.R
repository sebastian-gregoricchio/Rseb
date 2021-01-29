#' @title Volcano plot generator for RNA-seq data.
#'
#' @description Generates a volcano plot in order to visualize the differentially expressed genes. The plot is highly customizable.
#'
#' @param log2FC_data Numeric vector containing the log2(FoldChange) values of each gene.
#' @param padj_data Numeric vector of p-values. Use of adjusted p-values is recommended.
#' @param FC_t Value of the threshold to use for the fold change expression to define differentially expressed genes, expressed as linear value. By default 1.5 and by consequence 1/1.5.
#' @param p_t Value of the threshold to use for the p-values to define differentially expressed genes, expressed as linear value. By default 0.05.
#' @param FC_unresponsive_rigth Value of the threshold to use for the fold change expression to define unresponsive genes when \code{FC > 1}, expressed as linear value. By default 1.1. If \code{NULL} it will be calculated symmetrically from \code{FC_NoResp_left} as 1/\code{FC_NoResp_left}.
#' @param FC_unresponsive_left Value of the threshold to use for the fold change expression to define unresponsive genes when \code{FC < 1}, expressed as linear value. By default \code{1/FC_unresponsive_rigth}. If \code{NULL} it will be calculated symmetrically from \code{FC_NoResp_rigth} as 1/\code{FC_NoResp_rigth}.
#' @param x_ends Numeric positive value to define manually the range of the X-axis: it will be calculated as \code{c(-x_ends, x_ends)}, for this reason the plot will be symmetrical. By default \code{NULL}, the range is assigned automatically and the plot can be asymmetrical.
#' @param y_min Numeric value for the minimum value of the Y-axis. By default 0. Set it to \code{NULL} for automatic computation.
#' @param y_max Numeric value for the maximum value of the Y-axis. By default \code{NULL}.
#' @param left_label String to indicate the label to use for the set of genes in the left side of the graph (those with \code{FoldChange < 1/FC_t} and \code{p.value < p_t}. By default \code{"UP"}.
#' @param right_label String to indicate the label to use for the set of genes in the right side of the graph (those with \code{FoldChange > FC_t} and \code{p.value < p_t}. By default \code{"DOWN"}.
#' @param unresponsive_label String to indicate the label to use for the set of unresponsive genes (those with \code{FC_unresponsive_left < FoldChange < FC_unresponsive_rigth} and \code{p.value > p_t}. By default \code{"NoResp"}.
#' @param null_label String to indicate the label to use for the set of null genes (those with \code{1/FC_t < FoldChange < FC_t} and \code{p.value < p_t}. By default \code{"NULL"}.
#' @param names String vector with the names to be plotted if required, eg. gene names. By default \code{as.character(c(1:length(log2FC_data)))}.
#' @param left_names Logic value to indicate if to print the set of differentially expressed genes in the left side of the graph (those with \code{FoldChange < 1/FC_t} and \code{p.value < p_t}. By default \code{FALSE}.
#' @param right_names Logic value to indicate if to print the set of differentially expressed genes in the right side of the graph (those with \code{FoldChange > FC_t} and \code{p.value < p_t}. By default \code{FALSE}.
#' @param padding Logic value to indicate if to plot the padding around the names of genes. By default \code{FALSE}.
#' @param names_size Numeric value to define de size of the point names size. By default 10.
#' @param print_plot Logic value to define whether to print the volcano plot once created. By default \code{FALSE}.
#' @param left_color String to indicate the color to use for the set of genes in the left side of the graph (those with \code{FoldChange < 1/FC_t} and \code{p.value < p_t}. By default \code{"#00BA38"}, a green.
#' @param right_color String to indicate the color to use for the set of genes in the right side of the graph (those with \code{FoldChange > FC_t} and \code{p.value < p_t}. By default \code{"#F8766D"}, a pink/red.
#' @param unresponsive_color String to indicate the color to use for the set of unresponsive genes (those with \code{FC_unresponsive_left < FoldChange < FC_unresponsive_rigth} and \code{p.value > p_t}. By default \code{"#00A5CF"}, a light blue.
#' @param null_color String to indicate the color to use for the set of null genes (those with \code{1/FC_t < FoldChange < FC_t} and \code{p.value < p_t}. By default \code{"gray30"}, a dark gray.
#' @param point_size Numeric value to define de size of the points. By default 0.5.
#' @param legend Logic value to define if to print the legend. By default \code{TRUE}.
#' @param legend_title A string to indicate the label of the legend title. By default \code{"Expression status"}.
#' @param x_label A string to indicate the X-axis label. By default \code{"log2(fold change expression)"}.
#' @param y_label A string to indicate the Y-axis label. By default \code{"-log10(p-value adjusted)"}.
#' @param title A string to indicate the title of the plot. By default \code{"Volcano plot"}.
#' @param sub_title A string to indicate the subtitle of the plot. By default \code{NULL}, no subtitle is written.
#' @param add_threshold_lines Logic value to define if lines for the thresholds, both FC and p.value, should be plotted. By default \code{TRUE}.
#' @param threshold_line_color String to define the color of the threshold lines. By default \code{"gray70"}
#' @param threshold_line_type String or numeric value to define the threshold lines type. Both numeric and string standard R codes are accepted. By default \code{"dotted"}, equivalent to \code{2}.
#' @param font_family String to define the font family to use in the plot writings. By default \code{"Helvetica"}.
#' @param font_size Numeric value to define the font size. By default 12.
#'
#' @return A plot in ggplot2 format.
#'
#' @export volcano
#'
# @import ggplot2
# @import ggrepel

volcano = function(log2FC_data,
                        padj_data,

                        #set the thresholds
                        FC_t = 1.5,
                        p_t = 0.05,
                        FC_unresponsive_rigth = 1.1,
                        FC_unresponsive_left = 1/FC_unresponsive_rigth,

                        # axes values
                        x_ends = NULL, # x extremities: the plot will be symmetric
                        y_min = 0,
                        y_max = NULL,

                        # categorie names
                        left_label = "UP",
                        right_label = "DOWN",
                        unresponsive_label = "NoResp",
                        null_label = "NULL",

                        # genes names
                        names = as.character(c(1:length(log2FC_data))),
                        left_names = FALSE,
                        right_names = FALSE,
                        padding = FALSE,
                        names_size = 10,

                        # output option
                        print_plot = F,


                        #--------------------------------------------------
                        ## PLOT SETTINGS
                        # colors
                        left_color = "#00BA38",     # green
                        right_color = "#F8766D",   # red
                        unresponsive_color = "#00A5CF", # blue
                        null_color = "gray30",

                        # points size
                        point_size = 0.5,

                        # legend
                        legend = TRUE,
                        legend_title = "Expression status",

                        # lables
                        x_label = bquote('log'['2']*'(Fold Change expression)'),
                        y_label = bquote('-log'['10']*'(p-value'['adjusted']*')'),

                        # title
                        title = "Volcano plot",
                        sub_title = NULL,

                        # reference lines
                        add_threshold_lines = T,
                        threshold_line_color = "gray70",
                        threshold_line_type = "dotted",

                        # font parameters
                        font_family = "Helvetica",
                        font_size = 12) {

  if (length(log2FC_data) != length(padj_data)) {
    return(warning("ERROR: length(log2FC_data) is different from length(padj_data)."))
  }

  # Generate a table containing the data and their status
  table = data.frame(FC = log2FC_data,
                     padj = padj_data,
                     DE_status = factor(Rseb::DE.status(log2FC = log2FC_data,
                                                        p.value.adjusted = padj_data,
                                                        FC_threshold = FC_t,
                                                        FC_NoResp_left = FC_unresponsive_left,
                                                        FC_NoResp_rigth = FC_unresponsive_rigth,
                                                        p.value_threshold = p_t,
                                                        low.FC.status.label = left_label,
                                                        high.FC.status.label = right_label,
                                                        unresponsive.label = unresponsive_label,
                                                        null.label = null_label),
                                        levels = c(left_label, right_label, unresponsive_label, null_label)),
                     ID = names)

  # create the plot
  require(ggplot2)
  p =
    ggplot(table,
           aes(x = FC,
               y = -log10(padj),
               col = DE_status)) +
    geom_point(size = point_size) +
    ggtitle(label = title,
            subtitle = sub_title) +
    scale_color_manual(legend_title,
                       values = c(left_color,
                                  right_color,
                                  unresponsive_color,
                                  null_color)) +
    xlab(x_label) +
    ylab(y_label) +
    theme_classic() +
    theme(text = element_text(family = font_family, size = font_size),
          axis.text = element_text(color = "#000000"),
          axis.title = element_text(color = "#000000"),
          axis.line = element_line(color = "#000000"),
          axis.ticks = element_line(color = "#000000"),
          title = element_text(color = "#000000"),
          legend.title = element_text(color = "#000000"),
          legend.text = element_text(color = "#000000"))


  # remove legend if necessary
  if (legend == F) {p = p + theme(legend.position = "none")}

  # get axes limits
  if (is.null(x_ends)) {
    x.limit = ceiling(max(abs(ggplot_build(p)$layout$panel_params[[1]]$x$limits)))
  } else {x.limit = x_ends}
  build.plot = p +
    xlim(c(-x.limit, x.limit)) +
    ylim(c(y_min,
           ifelse(test = is.null(y_max),
                  yes = ceiling(max(abs(ggplot_build(p)$layout$panel_params[[1]]$y$limits))),
                  no = y_max)))

  build = ggplot_build(build.plot)

  # add thresholds lines and correspondent ticks
  if (add_threshold_lines == T) {
    p = p +

      # plot the p-value threshold
      geom_hline(yintercept = -log10(p_t),
                 linetype = threshold_line_type,
                 col = threshold_line_color) +

      # plot FC thresholds
      geom_vline(xintercept = c(-log2(FC_t), log2(FC_t)),
                 linetype = threshold_line_type,
                 col = threshold_line_color)


    # Add ticks of thresholds
    p = p +
      scale_y_continuous(breaks = sort(c(build$layout$panel_params[[1]]$y.sec$breaks,
                                         round(-log10(p_t), 1))),
                         limits = c(y_min,
                                    ifelse(test = is.null(y_max),
                                           yes = ceiling(max(abs(ggplot_build(p)$layout$panel_params[[1]]$y$limits))),
                                           no = y_max))) +
      scale_x_continuous(breaks = sort(c(build$layout$panel_params[[1]]$x.sec$breaks,
                                         round(log2(FC_t), 1), round(-log2(FC_t), 1))),
                         limits = c(-x.limit, x.limit))
  }

  # Show names of up/left genes
  if (left_names == T) {
    require(ggrepel)
    if (padding == T) {
    p = p +
      geom_label_repel(data = subset(table, table$DE_status == left_label),
                       aes(label = ID),
                       box.padding   = 0.35,
                       point.padding = 0.5,
                       segment.color = 'grey50',
                       size = names_size,
                       show.legend = F)}
    if (padding == F) {
      p = p +
        geom_text_repel(data = subset(table, table$DE_status == left_label),
                        aes(label = ID),
                        segment.color = 'grey50',
                        size = names_size,
                        show.legend = F)}
    }

  # Show names of down/right genes
  if (right_names == T) {
    require(ggrepel)
    if (padding == T) {
    p = p +
      geom_label_repel(data = subset(table, table$DE_status == right_label),
                       aes(label = ID),
                       box.padding   = 0.35,
                       point.padding = 0.5,
                       segment.color = 'grey50',
                       size = names_size,
                       show.legend = F)}
    if (padding == F) {
      p = p +
        geom_text_repel(data = subset(table, table$DE_status == right_label),
                        aes(label = ID),
                        segment.color = 'grey50',
                        size = names_size,
                        show.legend = F)}
    }


  # output and printing
  if (print_plot == T) {
    print(p)
    return(p)
  } else {return(p)}

} # END function
