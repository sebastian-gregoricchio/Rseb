#' @title qPCR RNA expression experimental replicates mean calculator.
#'
#' @description This function allows to generate a table and a plot result of the mean of different replicates of an experiment starting from analyses performed by \link{qPCR.rna.exp}.
#'
#' @param reps.list A list of \link{qPCR.rna.exp} results (and/or 'mean_FC_housekeeping' data.frames).
#' @param housekeeping.genes String vector with the list of genes that have to be used as target genes. By default \code{NULL}: if the input is a list of qPCR.rna.exp objects, the housekeeping genes are retrieved automatically.
#' @param exclude.samples String vector indicating the samples that should be exuded in the expression and FoldChange plots. By default \code{NULL}.
#' @param exclude.housekeeping.genes Logic value to indicate whether the housekeeping genes should be excluded in the plot. By default \code{TRUE}.
#' @param plot.color Single string to indicate the color to use for the bar plot. Default value: \code{#D1718B}.
#' @param fix.y.axis Logic value indicating whether the y-axis of the plots should be kept fixed among all the genes. By default \code{FALSE}.
#' @param text.size Numeric value to indicate the size of the text for the number above the bars. Default \code{3}.
#' @param x.labels.rotation Numeric value indicating the degrees of x-axis's labels rotation. By default \code{45}.
#' @param force Logic value to indicate whether the analyses should be performed also when the reference sample is not the same among the replicates. By default \code{FALSE}.

#'
#' @return The function returns a list containing:
#' \itemize{
#'   \item \code{mean.reps.table}: a data.frame containing the mean, number of reps (n), SD and SEM for each sample-target combination;
#'   \item \code{mean.reps.FC.plot}: a plot showing the replicates mean FoldChange expression over the reference Sample (facet_wrapped by gene).
#'  }
#'
#'
#' @export qPCR.rna.mean.reps



qPCR.rna.mean.reps = function(
    reps.list,
    housekeeping.genes = NULL,
    exclude.samples = NULL,
    exclude.housekeeping.genes = TRUE,
    plot.color = "#d1718b",
    fix.y.axis = FALSE,
    text.size = 3,
    x.labels.rotation = 45,
    force = F) {


  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # libraries
  require(dplyr)
  require(ggplot2)


  # Check plot color
  if ((plot.color != "#d1718b") & !Rseb::is.color(plot.color)) {
    message("The plot.color value is not recognized. Coerced to '#d1718b' by default.")
    plot.color = "#d1718b"
  }


  # Select columns in the list
  check_names = c()
  housekeeping.genes.list = c()

  for (i in 1:length(reps.list)) {
    if ("list" %in% class(reps.list[[i]])) {
      if (is.null(housekeeping.genes)) {housekeeping.genes.list = c(housekeeping.genes.list, names(reps.list[[i]]$expression.plots))}
      reps.list[[i]] = reps.list[[i]]$analyzed.data$mean_FC_housekeeping
    }
    reps.list[[i]] = reps.list[[i]][, c(1,2,(ncol(reps.list[[i]])-1))]
    check_names = c(check_names, colnames(reps.list[[i]][3]))
  }

  housekeeping.genes.list = unique(housekeeping.genes.list)


  # Check excluded samples
  if (!is.null(exclude.samples)) {
    if (F %in% (exclude.samples %in% reps.list[[i]]$`Sample Name`)) {
      message(paste0("The 'exclude.samples' values are not among the sample_IDs present in your file -> ",
                     paste(unique(reps.list[[i]]$`Sample Name`), collapse = ", "), ". \nHowever, the function will run anyway."))
    }
  }


  # Define order
  samples_order = unique(reps.list[[i]]$`Sample Name`)


  # Check that the FC are coherent
  if ((length(unique(check_names)) > 1) & (force == F)) {
    return(warning("The FoldChanges provided correspond to comparisons made on different reference samples.\nIf you wish to perform the analyses anyway, set the parameter 'force = TRUE'."))
  }


  # Combine tables
  combined_table = do.call("rbind", reps.list)
  comparison_name = colnames(combined_table)[3]
  colnames(combined_table)[3] = "FoldChange"

  # Calculate the statistics
  stat_table =
    combined_table %>%
    dplyr::group_by(`Sample Name`, `Target Name`) %>%
    dplyr::summarise(n = n(),
                     mean_FC = mean(FoldChange, na.rm = T),
                     SD = sd(FoldChange, na.rm = T),
                     .groups = "drop_last") %>%
    dplyr::mutate(SEM = SD / sqrt(n)) %>%
    dplyr::filter(!(`Sample Name` %in% exclude.samples))

  stat_table =
    stat_table %>%
    dplyr::mutate(`Sample Name` = factor(`Sample Name`, samples_order)) %>%
    dplyr::arrange(`Target Name`, `Sample Name`)


  # Filter housekeeping if required
  if (exclude.housekeeping.genes == TRUE) {
    table = stat_table %>% dplyr::filter(!(`Target Name` %in% housekeeping.genes.list))
  } else (table = stat_table)

  table =
    dplyr::mutate(.data = table, `Sample Name` = factor(`Sample Name`, samples_order)) %>%
    dplyr::arrange(`Target Name`, `Sample Name`) %>%
    dplyr::filter(!(`Sample Name` %in% exclude.samples))


  # Generate plot
  mean_FC_plot =
    ggplot(data = table,
           aes(x = `Sample Name`,
               y = mean_FC)) +
    geom_bar(stat = "identity",
             fill = plot.color,
             position = position_dodge(width=0.9)) +
    geom_errorbar(aes(ymin = mean_FC-SEM,
                      ymax = mean_FC+SEM),
                  position = position_dodge(width=0.9),
                  width = 0.3,
                  color = "#000000") +
    geom_text(data = table,
              aes(x = `Sample Name`,
                  y = mean_FC+SEM,
                  label = round(mean_FC,2)),
              position = position_dodge(width=0.9),
              vjust=-0.25,
              size = text.size,
              inherit.aes = F) +
    geom_text(data = table,
              aes(x = `Sample Name`,
                  y = 0,
                  label = n),
              position = position_dodge(width=0.9),
              vjust=-0.25,
              size = text.size-0.5,
              inherit.aes = F) +
    facet_wrap(~`Target Name`,
               scales = ifelse(test = fix.y.axis == T,
                               yes = "fixed",
                               no = "free")) +
    ylab("mean rep FoldChange \U00B1 SEM") +
    ggtitle("mean FoldChange of all replicates") +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = 2, color = "#000000") +
    theme(axis.text = element_text(color = "#000000"),
          axis.text.x = element_text(angle = x.labels.rotation,
                                     hjust = 1))



  # Return output
  return(list(mean.reps.table = stat_table,
              mean.reps.FC.plot = mean_FC_plot))

} # END function


