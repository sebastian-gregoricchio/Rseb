#' @title qPCR RNA expression experimental replicates mean calculator.
#'
#' @description This function allows to generate a table and a plot result (FoldChange and normalized Expression) of the mean of different replicates of an experiment starting from analyses performed by \link{qPCR.rna.exp}.
#'
#' @param reps.list A list of \link{qPCR.rna.exp} result objects.
#' @param reference.sample Single string indicating the name of the sample to use as reference for the computation of the FoldChanges. By default \code{NULL}: if the input is a list of qPCR.rna.exp objects, the reference.sample is retrieved automatically. However, if the number of reference used are multiple and/or not shared among replicates, the first sample in the order is used as reference.
#' @param exclude.samples String vector indicating the samples that should be exuded in the expression and FoldChange plots. By default \code{NULL}.
#' @param exclude.housekeeping.genes Logic value to indicate whether the housekeeping genes should be excluded in the plot. By default \code{TRUE}.
#' @param plot.color Single string to indicate the color to use for the bar plot. Default value: \code{#D1718B}.
#' @param fix.y.axis Logic value indicating whether the y-axis of the plots should be kept fixed among all the genes. By default \code{FALSE}.
#' @param text.size Numeric value to indicate the size of the text for the number above the bars. Default \code{3}.
#' @param x.labels.rotation Numeric value indicating the degrees of x-axis's labels rotation. By default \code{45}.

#'
#' @return The function returns a list containing:
#' \itemize{
#'   \item \code{mean.reps.data.table}: a list of data.frames, one per housekeeping gene and the mean of all housekeeping genes, containing the number of reps (n), SD and SEM for each sample-target combination for both normalized expression and FoldChanges;
#'   \item \code{mean.reps.exp.plots}: a list of a plots, one per housekeeping gene, showing the replicates' normalized expression over the reference Sample (facet_wrapped by gene);
#'   \item \code{mean.reps.FC.plots}: a list of a plots, one per housekeeping gene and the mean of all housekeeping genes, showing the replicates' mean FoldChange expression over the reference Sample (facet_wrapped by gene).
#'  }
#'
#'
#' @export qPCR.rna.mean.reps



qPCR.rna.mean.reps = function(
  reps.list,
  reference.sample = NULL,
  exclude.samples = NULL,
  exclude.housekeeping.genes = TRUE,
  plot.color = "#d1718b",
  fix.y.axis = FALSE,
  text.size = 3,
  x.labels.rotation = 45) {


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
  housekeeping.genes.list = c()
  samples.list = c() # unique values
  all.samples = list() # a list with the unique ones for each table, need to find the common ones
  reference.sample.detected = c()

  for (i in 1:length(reps.list)) {
    if ("list" %in% class(reps.list[[i]])) {
      housekeeping.genes.list = unique(c(housekeeping.genes.list, names(reps.list[[i]]$expression.plots)))
      samples.list = unique(c(samples.list, unique(reps.list[[i]]$analyzed.data$mean_FC_housekeeping$`Sample Name`)))
      all.samples[[i]] = unique(reps.list[[i]]$analyzed.data$mean_FC_housekeeping$`Sample Name`)
      reference.sample.detected = unique(c(reference.sample.detected,
                                           gsub(pattern = "mean_FC_over_",
                                                replacement = "",
                                                x = names(reps.list[[1]]$mean_FC_housekeeping)[grep(pattern = "mean_FC_over_",
                                                                                                    x = names(reps.list[[1]]$mean_FC_housekeeping))])))
      reps.list[[i]] = reps.list[[i]]$analyzed.data
    }
  }


  # Check reference sample
  if (!is.null(reference.sample) | (length(reference.sample) > 1)) {
    if (!(reference.sample %in% Reduce(intersect, all.samples))) {
      return(warning(paste0("The 'reference.sample' must be among the sample_IDs present in ALL your files -> ",
                            paste(unique(results$`Sample Name`), collapse = ", "))))
    }
  } else if (length(reference.sample.detected) == 1) {
    reference.sample = reference.sample.detected[1]
    message(paste0("The 'reference.sample' has been automatically detected and set to: '", reference.sample, "'."))
  } else {
    reference.sample = Reduce(intersect, all.samples)[1]
    message(paste0("The 'reference.sample' has been automatically set to: '", reference.sample, "'."))
  }


  # Check excluded samples
  if (!is.null(exclude.samples)) {
    if (F %in% (exclude.samples %in% samples.list)) {
      message(paste0("The 'exclude.samples' values are not among the sample_IDs present in your file -> ",
                     paste(unique(samples.list), collapse = ", "), ". \nHowever, the function will run anyway."))
    }
  }


  # Define the sample order
  samples.order = c(reference.sample, samples.list[samples.list != reference.sample])


  # --------------------------------------------------------------------------------------
  # Combine the tables by housekeeping gene
  housekeeping.table.list = list()

  for (i in 1:length(housekeeping.genes.list)) {
    housekeeping.table.list[[i]] = list()

    for (k in 1:length(reps.list)) {
      if (housekeeping.genes.list[i] %in% names(reps.list[[k]])) {
        current_tb = reps.list[[k]][[housekeeping.genes.list[i]]][grepl(pattern = "(Sample|Target Name)|exp_norm",
                                                                        x = names(reps.list[[k]][[housekeeping.genes.list[i]]]))]

        # Recompute the FoldChange over the new reference.sample
        new_FC = list()
        for (z in 1:length(unique(current_tb$`Target Name`))) {
          tb = dplyr::filter(current_tb, `Target Name` == unique(current_tb$`Target Name`)[z])
          names(tb) = gsub(pattern = paste0("_to_",housekeeping.genes.list[i]), replacement = "", x = names(tb))
          new_FC[[z]] = dplyr::mutate(tb, FoldChange = exp_norm / tb$exp_norm[which(tb$`Sample Name` == reference.sample)])
        }
        new_FC = do.call(rbind, new_FC)
        housekeeping.table.list[[i]][[(length(housekeeping.table.list[[i]])+1)]] = new_FC
      }
    }

    # Merge reps table for the current housekeeping
    housekeeping.table.list[[i]] = do.call(rbind, housekeeping.table.list[[i]])

    # Compute stats for current housekeeping
    housekeeping.table.list[[i]] =
      housekeeping.table.list[[i]] %>%
      dplyr::filter(!is.na(exp_norm)) %>%
      dplyr::group_by(`Sample Name`, `Target Name`) %>%
      dplyr::summarise(n = n(),
                       mean_exp = mean(exp_norm, na.rm = T),
                       SD_exp = sd(exp_norm, na.rm = T),
                       mean_FC = mean(FoldChange, na.rm = T),
                       SD_FC = sd(FoldChange, na.rm = T),
                       .groups = "drop_last") %>%
      dplyr::mutate(SEM_exp = SD_exp / sqrt(n),
                    SEM_FC = SD_FC / sqrt(n)) %>%
      dplyr::filter(!(`Sample Name` %in% exclude.samples)) %>%
      dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = samples.order)) %>%
      dplyr::arrange(`Target Name`, `Sample Name`)

    housekeeping.table.list[[i]] = Rseb::move.df.col(data.frame = housekeeping.table.list[[i]],
                                                     move.command = "SEM_exp after SD_exp")

  }
  names(housekeeping.table.list) = housekeeping.genes.list


  # --------------------------------------------------------------------------------------
  # Combine the tables for FoldChanges of all housekeeping together
  FC_all_housekeeping = list()

  for (i in 1:length(reps.list)) {
    FC_current_rep = list()

    # Calculating FoldChange for each housekeeping before to make the mean all together
    for (k in 1:(length(reps.list[[i]])-1)) {
      temp_housekeeping = data.frame(reps.list[[i]][[k]], check.names = F)[grepl(pattern = "(Sample|Target Name)|exp_norm", x = names(reps.list[[i]][[k]]))]

      # Recompute the FoldChange over the new reference.sample for the current housekeeping
      new_FC = list()
      for (z in 1:length(unique(temp_housekeeping$`Target Name`))) {
        tb = dplyr::filter(temp_housekeeping, `Target Name` == unique(temp_housekeeping$`Target Name`)[z])
        names(tb) = gsub(pattern = "_to_.*$", replacement = "", x = names(tb))
        new_FC[[z]] = dplyr::mutate(tb, FoldChange = exp_norm / tb$exp_norm[which(tb$`Sample Name` == reference.sample)])
      }
      new_FC = do.call(rbind, new_FC)
      FC_current_rep[[k]] = new_FC
    }

    FC_all_housekeeping[[i]] = do.call(rbind, FC_current_rep)

    FC_all_housekeeping[[i]] =
      FC_all_housekeeping[[i]] %>%
      dplyr::group_by(`Sample Name`, `Target Name`) %>%
      dplyr::summarise(mean_FoldChange = mean(FoldChange, na.rm = T),
                       .groups = "drop_last")
  }

  FC_all_housekeeping = do.call(rbind, FC_all_housekeeping)

  FC_all_housekeeping =
    FC_all_housekeeping %>%
    dplyr::filter(!is.na(mean_FoldChange)) %>%
    dplyr::group_by(`Sample Name`, `Target Name`) %>%
    dplyr::summarise(n = n(),
                     mean_FC = mean(mean_FoldChange, na.rm = T),
                     SD_FC = sd(mean_FoldChange, na.rm = T),
                     .groups = "drop_last") %>%
    dplyr::mutate(SEM_FC = SD_FC / sqrt(n)) %>%
    dplyr::filter(!(`Sample Name` %in% exclude.samples)) %>%
    dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = samples.order)) %>%
    dplyr::arrange(`Target Name`, `Sample Name`)


  # --------------------------------------------------------------------------------------
  # Generate a final results table
  all_results_table = Rseb::combine.lists(list(housekeeping.table.list, list(FC_all_housekeeping)))
  names(all_results_table)[length(all_results_table)] = "mean_FC.means"


  # --------------------------------------------------------------------------------------
  # Generate expression plot
  expression_plots = list()

  for (i in 1:length(housekeeping.table.list)){

    # Filter housekeeping if required
    if (exclude.housekeeping.genes == TRUE) {
      table = housekeeping.table.list[[i]] %>% dplyr::filter(!(`Target Name` %in% housekeeping.genes.list))
    } else (table = housekeeping.table.list[[i]])

    # Build the plot expression plot
    expression_plots[[i]] =
      ggplot(data = table,
             aes(x = `Sample Name`,
                 y = mean_exp)) +
      geom_bar(stat = "identity",
               fill = plot.color,
               position = position_dodge(width=0.9)) +
      geom_errorbar(aes(ymin = mean_exp-SEM_exp,
                        ymax = mean_exp+SEM_exp),
                    position = position_dodge(width=0.9),
                    width = 0.3,
                    color = "#000000") +
      geom_text(data = table,
                aes(x = `Sample Name`,
                    y = mean_exp+SEM_exp,
                    label = round(mean_exp,2)),
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
      ylab("mean normalized expression \U00B1 SEM") +
      ggtitle(paste("Exp. normalized to", names(housekeeping.table.list)[i],"- all reps")) +
      theme_classic() +
      theme(axis.text = element_text(color = "#000000"),
            axis.text.x = element_text(angle = x.labels.rotation,
                                       hjust = 1))

  }

  names(expression_plots) = names(housekeeping.table.list)


  # --------------------------------------------------------------------------------------
  # Generate FoldChange plots
  FoldChange_plots = list()

  for (i in 1:length(all_results_table)){

    # Filter housekeeping if required
    if (exclude.housekeeping.genes == TRUE) {
      table = all_results_table[[i]] %>% dplyr::filter(!(`Target Name` %in% housekeeping.genes.list))
    } else (table = all_results_table[[i]])

    # Build the FC expression plot
    FoldChange_plots[[i]] =
      ggplot(data = table,
             aes(x = `Sample Name`,
                 y = mean_FC)) +
      geom_bar(stat = "identity",
               fill = colorspace::darken(plot.color),
               position = position_dodge(width=0.9)) +
      geom_errorbar(aes(ymin = mean_FC-SEM_FC,
                        ymax = mean_FC+SEM_FC),
                    position = position_dodge(width=0.9),
                    width = 0.3,
                    color = "#000000") +
      geom_text(data = table,
                aes(x = `Sample Name`,
                    y = mean_FC+SEM_FC,
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
      ylab("mean FoldChange nor. exp. \U00B1 SEM") +
      ggtitle(paste0("mean FC norm. exp. over '", reference.sample, "' (", names(all_results_table)[i],") - all reps")) +
      theme_classic() +
      geom_hline(yintercept = 1, linetype = 2, color = "#000000") +
      theme(axis.text = element_text(color = "#000000"),
            axis.text.x = element_text(angle = x.labels.rotation,
                                       hjust = 1))

  }

  names(FoldChange_plots) = names(all_results_table)


  # --------------------------------------------------------------------------------------
  # Return output
  return(list(mean.reps.data.table = all_results_table,
              mean.reps.exp.plots = expression_plots,
              mean.reps.FC.plots = FoldChange_plots))

} # END function
