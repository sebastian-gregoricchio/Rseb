#' @title Incucyte analysis tool
#'
#' @description This tool generates confluency plots over time starting from incucyte data.
#'
#' @param metadata Path to a table or a data.frame with at least two columns: 'well.ID' and 'group'. The 'well.id' indicates the ID of the wells/columns in the rawData table; the 'group' indicates the name of the group at which each column belongs (wells belonging to the same column will be averaged together). A third column 'color' can be used to indicate the color; notice that each group should have the same color.
#' @param raw.data Path to the rawData table or a data.frame. First column must be 'Elapsed' (the time in hours) and the other columns the wells (e.g., A1, A2, B5, B8, ...).
#' @param comparisons List of vectors containing the groups to be used in each comparison: one vector per comparison. E.g.: list(c("group_A", "group_B"), c("group_B", "group_C", "group_E")). Default: \code{"all"}.
#' @param error String to indicate the error type. Possible choices: 'SEM', 'SD'. Default: \code{SEM}.
#' @param normalization.method String to indicate the normalization method. Possible choices: "none", "division", "subtraction". Default: \code{none}. Division: for each group the values are divided by the first timepoint.
#' @param start.hour Number of hours to be excluded from the calculations. Default: \code{0}.
#' @param plot.days Logical value to indicate whether to plotted days instead of hours. Default: \code{FALSE}.
#' @param show.error.ribbon Logical value to indicate whether to show the error ribbon. Default: \code{TRUE}.
#' @param show.error.bars Logical value to indicate whether to show the error bars. Default: \code{FALSE}.
#' @param show.points Logical value to indicate whether to show the individual average points at each time point. Default: \code{TRUE}.
#' @param show.lines Logical value to indicate whether to show the interpolation line. Default: \code{TRUE}.
#' @param show.legend Logical value to indicate whether to show the color legend. Default: \code{TRUE}.
#' @param same.y.scale Logical value to indicate whether the y-axis should show the same limit range. Default: \code{TRUE}.
#' @param group.order String vector with the specific order to use for the groups indicated in the metadata (factor levels). Default: \code{NULL}.
#' @param error.transparency Number between 0 and 1 indicating the transparency (alpha) to use for the error ribbon. Default: \code{0.25}.
#' @param point.size Number for the point size. Default: \code{1}.
#' @param line.type A vector indicating the line type to use for each group (both numeric and string values accepted). Default: \code{1} (applied to all the groups).
#' @param line.smooth.span Numeric value of the 'span' for the smoothing of the interpolation line. Default: \code{0.25}.
#' @param colors A vector indicating the color to use for each group (any R color format is accepted). Default: \code{NULL} (random colors are generated).
#' @param skip.head.lines.in.data Number of lines to be skipped at the beginning of the raw.data file. Default: \code{1}.
#'
#'
#' @return The function returns a list containing:
#' \itemize{
#'   \item \code{metadata}: table used for the metadata;
#'   \item \code{raw.data}: tables used as raw.data;
#'   \item \code{analyzed.data}: data.frame with the analyzed data (means, n, SEM, SD, groups, ...);
#'   \item \code{normalized.data}: data.frame with the normalized data used for the plotting;
#'   \item \code{plot.list}: names list of plots, one plot for each group plus 'all' comparison;
#'   \item \code{multiplot}: a plot with all the plots in the plot.list .
#'  }
#'
#' @export incucyte



incucyte =
  function(metadata,
           raw.data,
           comparisons = "all",
           error = "SEM",
           normalization.method = "none",
           start.hour = 0,
           plot.days = FALSE,
           show.error.ribbon = TRUE,
           show.error.bars = FALSE,
           show.points = TRUE,
           show.lines = TRUE,
           show.legend = TRUE,
           same.y.scale = TRUE,
           group.order = NULL,
           error.transparency = 0.25,
           point.size = 1,
           line.type = 1,
           line.smooth.span = 0.25,
           colors = NULL,
           skip.head.lines.in.data = 1) {

    #-----------------------------#
    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)
    #-----------------------------#

    # Libraries
    require(dplyr)
    require(ggplot2)


    # Check normalization method
    norm.methods = c("none", "division", "subtraction")
    if (!(tolower(normalization.method) %in% norm.methods)) {
      return(warning(paste0("Normalization method not supported. Possibile choices: '",
                            paste(norm.methods, collapse = "', '"), "'.")))
    }


    # Hours or days?
    time.factor = ifelse(test = isTRUE(plot.days), yes = 24, no = 1)
    time.type = ifelse(test = isTRUE(plot.days), yes = "days", no = "hours")


    # reading metadata
    if ("character" %in% class(metadata)) {
      metadata = data.frame(data.table::fread(metadata, blank.lines.skip = T))
    } else if (!("data.frame" %in% class(metadata))) {
      return(warning("The metadata table must be a data.frame or a path to a table."))
    }

    if (ncol(metadata) > 3) {metadata = metadata[,1:3]}
    colnames(metadata) = c("well.ID", "group", "color")[1:ncol(metadata)]



    # Define order
    if (is.null(group.order)) {
      group.order = unique(metadata$group)
    } else {
      if (length(unique((unique(metadata$group) %in% group.order))) > 1) {
        return(warning(paste0("In the group.order there groups/samples missing if compared to the metadata provided.")))
      }
    }



    # Define colors
    if (is.null(colors) & ncol(metadata)<3) {
      n.colors = length(unique(metadata$group))
      color.list = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      random.colors = color.list[runif(n = n.colors, min = 1, max = length(color.list))]
      color.vector = random.colors
      names(color.vector) = unique(metadata$group)
    } else if (is.null(colors) & ncol(metadata) == 3){
      color.tb = unique(metadata[,2:3])
      color.vector = color.tb[,3]
      names(color.vector) = color.tb[,2]
    } else if (!is.null(colors)) {
      if (ncol(metadata) == 3) {
        message("Color list is provided but other colors are defined in the metadata table. Only the first set will be used.")
      }
      if (is.null(names(colors))) {
        names(colors) = group.order
      } else if (length(unique(unique(metadata[,2]) %in% names(colors))) > 1) {
        return(warning("The names provided in the color list do not include all the samples/groups indicated.\n"))
      }
      color.vector = colors
    }



    # Define line.type
    if (length(line.type) > 1) {
      if (is.null(names(line.type))) {
        names(line.type) = group.order
      } else if (length(unique(unique(metadata[,2]) %in% names(line.type))) > 1) {
        return(warning("The names provided in the line.type list do not include all the samples/groups indicated.\n"))
      }
    } else if (length(line.type) == 1) {
      line.type = rep(line.type, length(group.order))
      names(line.type) = group.order
    }



    # Check comparisons
    if (comparisons[1] != "all") {
      if ("list" %in% class(comparisons)) {
        for (i in 1:length(comparisons)) {
          if (length(unique((comparisons[[i]] %in% unique(metadata$group)))) > 1) {
            return(warning(paste0("The element #",i,"of the comparison list contains comparisons between groups not present among the groups provided.")))
          }
          comparisons[[i]] = unique(comparisons[[i]])
          names(comparisons)[i] = paste(comparisons[[i]], collapse="_")
        }
        comparisons = Rseb::combine.lists(list(list(unique(metadata$group)),
                                               comparisons))
        names(comparisons)[1] = "all.groups"
      } else {
        for (i in 1:length(comparisons)) {
          if (length(unique((comparisons[i] %in% unique(metadata$group)))) > 1) {
            return(warning(paste0("The element #",i,"of the comparison list contains comparisons between groups not present among the groups provided.")))
          }
        }
        comparisons = list(unique(metadata$group), comparisons)
        names(comparisons) = c("all.groups", paste(comparisons[[2]], collapse="_"))
      }
    }



    # Importing data
    if ("character" %in% class(raw.data)) {
      raw.data = data.table::fread(raw.data, blank.lines.skip = T, skip = skip.head.lines.in.data)
    } else if (!("data.table" %in% class(metadata))) {
      return(warning("The metadata table must be a data.frame or a path to a table."))
    }

    raw.data = dplyr::mutate(raw.data,
                             across(1:ncol(raw.data),
                                    function(x){gsub(pattern = ",",
                                                     replacement = ".",
                                                     x = x)}))

    raw.data = dplyr::mutate(raw.data, across(2:ncol(raw.data), as.numeric))



    # Compute the mean by group
    groups = unique(metadata$group)

    analyzed.tb = data.frame()
    for (i in 1:length(groups)) {
      ids = (dplyr::filter(metadata, group %in% groups[i]))$well.ID
      subset.tb = as.data.frame(raw.data)[,c(1:2,grep(paste(ids,collapse = "|"), colnames(raw.data)))]
      subset.tb =
        subset.tb %>%
        dplyr::mutate(group = groups[i],
                      n = ncol(subset.tb)-2,
                      mean = rowMeans(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T),
                      SD = matrixStats::rowSds(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T)) %>%
        dplyr::mutate(SEM = SD / sqrt(n))

      analyzed.tb = rbind(analyzed.tb,
                          subset.tb[,c(1:2, (ncol(subset.tb)-4):ncol(subset.tb))])
    }


    ######################## NORMALIZATION ##############################
    normalized.tab.all = data.frame()

    #NONE
    if (tolower(normalization.method) == "none") {
      normalized.tab.all = analyzed.tb

    #DIVISION.FIRST.VALUE
    } else if (tolower(normalization.method) == "division") {
      raw.data.divided =
        dplyr::filter(raw.data, Elapsed >= start.hour) %>%
        dplyr::arrange(Elapsed)

      # division by the minimum value of each column
      first.value = unlist(raw.data.divided[1,3:ncol(raw.data.divided)])
      raw.data.divided[,3:ncol(raw.data.divided)] = data.frame(sweep(as.matrix(raw.data.divided[,3:ncol(raw.data.divided)]),
                                                                     2, first.value, FUN="/"))

      for (i in 1:length(groups)) {
        ids = (dplyr::filter(metadata, group %in% groups[i]))$well.ID
        subset.tb = as.data.frame(raw.data.divided)[,c(1:2,grep(paste(ids,collapse = "|"), colnames(raw.data.divided)))]
        subset.tb =
          subset.tb %>%
          dplyr::mutate(group = groups[i],
                        n = ncol(subset.tb)-2,
                        mean = rowMeans(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T),
                        SD = matrixStats::rowSds(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T)) %>%
          dplyr::mutate(SEM = SD / sqrt(n))

        normalized.tab.all = rbind(normalized.tab.all,
                                   subset.tb[,c(1:2, (ncol(subset.tb)-4):ncol(subset.tb))])
      }

    #SUBSTRACTION
    } else if (tolower(normalization.method) == "subtraction") {
      raw.data.subtraction =
        dplyr::filter(raw.data, Elapsed >= start.hour) %>%
        dplyr::arrange(Elapsed)

      # division by the minimum value of each column
      first.value = unlist(raw.data.subtraction[1,3:ncol(raw.data.subtraction)])
      raw.data.subtraction[,3:ncol(raw.data.subtraction)] = data.frame(sweep(as.matrix(raw.data.subtraction[,3:ncol(raw.data.subtraction)])+1,
                                                                             2, first.value, FUN="-"))

      for (i in 1:length(groups)) {
        ids = (dplyr::filter(metadata, group %in% groups[i]))$well.ID
        subset.tb = as.data.frame(raw.data.subtraction)[,c(1:2,grep(paste(ids,collapse = "|"), colnames(raw.data.subtraction)))]
        subset.tb =
          subset.tb %>%
          dplyr::mutate(group = groups[i],
                        n = ncol(subset.tb)-2,
                        mean = rowMeans(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T),
                        SD = matrixStats::rowSds(as.matrix(subset.tb[3:ncol(subset.tb)]), na.rm = T)) %>%
          dplyr::mutate(SEM = SD / sqrt(n))

        normalized.tab.all = rbind(normalized.tab.all,
                                   subset.tb[,c(1:2, (ncol(subset.tb)-4):ncol(subset.tb))])
      }
    }

    #####################################################################



    # select the error to plot
    if (toupper(error) == "SEM") {
      normalized.tab.all = dplyr::select(normalized.tab.all, -SD)
    } else {
      normalized.tab.all = dplyr::select(normalized.tab.all, -SEM)
    }
    colnames(normalized.tab.all)[6] = "error"



    # Generate plots
    comparison.plots = list()

    for (i in 1:length(comparisons)) {

      groups.to.use = if ("all" %in%comparisons[[i]]) {unique(normalized.tab.all$group)} else {comparisons[[i]]}

      normalized.tb = dplyr::filter(normalized.tab.all, group %in% groups.to.use)

      comparison.plots[[i]] =
        ggplot(data =
                 normalized.tb %>%
                 dplyr::mutate(group = factor(group, levels = group.order)) %>%
                 dplyr::filter(Elapsed >= start.hour),
               aes(x = Elapsed / time.factor,
                   y = mean,
                   color = group,
                   fill = group))

      if (isTRUE(show.error.ribbon)) {
        comparison.plots[[i]] =
          comparison.plots[[i]] +
          geom_ribbon(aes(ymin = mean - error,
                          ymax = mean + error),
                      alpha = error.transparency,
                      color = NA,
                      show.legend = show.legend) +
          scale_fill_manual(values = color.vector)
      }

      if (isTRUE(show.points)) {
        comparison.plots[[i]] =
          comparison.plots[[i]] +
          geom_point(size = point.size,
                     show.legend = show.legend) +
          scale_color_manual(values = color.vector)
      }

      if (isTRUE(show.error.bars)) {
        comparison.plots[[i]] =
          comparison.plots[[i]] +
          geom_errorbar(aes(ymin = mean - error,
                            ymax = mean + error),
                        show.legend = show.legend) +
          scale_color_manual(values = color.vector)
      }

      if (isTRUE(show.lines)) {
        comparison.plots[[i]] =
          comparison.plots[[i]] +
          geom_smooth(aes(linetype = group),
                      method = "loess",
                      formula = y ~ x,
                      span = line.smooth.span,
                      se = F,
                      show.legend = show.legend) +
          scale_linetype_manual(values = line.type) +
          scale_color_manual(values = color.vector)
      }


      comparison.plots[[i]] =
        comparison.plots[[i]] +
        ylab(bquote(atop(paste("Confluency (%) \u00B1 ", .(toupper(error))),
                         atop(paste("Normalization method: ", .(normalization.method)))))) +
        #ylab(paste0("Confluency (%) \u00B1 ", toupper(error))) +
        xlab(paste0("Time (", time.type ,")")) +
        ggtitle(paste(groups.to.use, collapse = " | ")) +
        theme_classic() +
        theme(axis.ticks = element_line(color = "#000000"),
              axis.text = element_text(color = "#000000"))
    }


    # Normalize y-axis if required
    plot.list = comparison.plots
    names(plot.list) = names(comparisons)

    if (isTRUE(same.y.scale) & length(plot.list)>1) {
      plot.list = Rseb::uniform.y.axis(plot.list)
    }



    # Export results
    colnames(normalized.tab.all)[6] = toupper(error)

    return(list(metadata = dplyr::left_join(metadata,
                                            dplyr::mutate(data.frame(color.vector, stringsAsFactors = F),
                                                          group = groups),
                                            by = "group"),
                raw.data = raw.data,
                analyzed.data = analyzed.tb,
                normalized.data = normalized.tab.all,
                plot.list = plot.list,
                multiplot = cowplot::plot_grid(plotlist = plot.list, align = "hv")))

  } # END function
