#' @title qPCR RNA expression analyses tool.
#'
#' @description Allows to easily analyse qPCR RNA expression data, including: technical replicates verification, gene expression normalization to housekeeping genes and FoldChanges over reference sample computation.
#'
#' @param results.file String indicating the full path to the results excel file or a data.frame containing at least the following columns: 'Sample Name', 'Target Name', 'CT'.
#' @param housekeeping.genes String vector with the list of genes that have to be used as target genes. By default \code{NULL}: an error message is printed.
#' @param max.delta.reps Numeric value indicating the maximum difference among replicate Ct. Default value: \code{0.5}.
#' @param reference.sample Single string indicating the name of the sample to use as reference for the computation of the FoldChanges. By default \code{NULL}: the first sample in the order is used as reference.
#' @param exclude.housekeeping.FC Logic value to indicate whether the housekeeping genes should be excluded in the FoldChanges plots. By default \code{TRUE}.
#' @param exclude.samples String vector indicating the samples that should be exuded in the expression and FoldChange plots. By default \code{NULL}.
#' @param fix.y.axis Logic value indicating whether the y-axis of the plots should be kept fixed among all the genes. By default \code{FALSE}.
#' @param x.labels.rotation Numeric value indicating the degrees of x-axis's labels rotation. By default \code{45}.
#' @param text.size Numeric value to indicate the size of the text for the number above the bars. Default \code{3}.
#' @param results.sheet.position Numeric value indicating the position of the results sheet in the excel file. by default \code{3}.
#' @param rows.to.skip How many rows must be skipped before to read the excel file. By default \code{44}.
#' @param file.header Logic value to indicate whether the results excel file contains an header. By default \code{TRUE}.
#' @param file.tail Logic value to indicate whether the results excel file contains extra rows at the end of the results. By default \code{TRUE}.
#' @param samples.order A string vector indicating all the samples in order. This order will be used to order the samples in the plots. By default \code{NULL}: the reference sample will be the first, the other will be kept in the order available in the results table.
#' @param ignore.reps.errors Logic value to define whether the difference between the Ct in replicates should be ignored: all the values are kept.
#'
#' @return The function returns a list containing:
#' \itemize{
#'   \item \code{original.table}: a data.frame containing the original results table;
#'   \item \code{reshaped.table}: a data.frame with the original results reorganized for the analyses;
#'   \item \code{reshaped.table.cleaned}: the reshaped data.frame upon filtering of the CT values (if required);
#'   \item \code{reps.validation.plot}: a plot representing a table with the differences two-by-two of the technical replicates (facet_wrapped by gene) where the cells have a red background if the difference is greater than the `max.delta.reps` value;
#'   \item \code{analyzed.data}: a named list of data.frames, one for each housekeeping gene and one for the foldChange mean of all housekeeping normalization, containing the normalized expression scores and the FoldChanges over the reference sample;
#'   \item \code{expression.plots}: a named list of plots, one for each housekeeping gene, showing the gene expression histograms (facet_wrapped by gene);
#'   \item \code{foldChange.plots}: a named list of plots, one for each housekeeping gene and one for the foldChange mean of all housekeeping normalization, showing the FoldChange expression over the reference Sample (facet_wrapped by gene).
#'  }
#'
#'
#' @export qPCR.rna.exp



qPCR.rna.exp = function(results.file,
                        housekeeping.genes = NULL,
                        max.delta.reps = 0.5,
                        reference.sample = NULL,
                        exclude.housekeeping.FC = TRUE,
                        exclude.samples = NULL,
                        fix.y.axis = FALSE,
                        x.labels.rotation = 45,
                        text.size = 3,
                        results.sheet.position = 3,
                        rows.to.skip = 44,
                        file.header = TRUE,
                        file.tail = TRUE,
                        samples.order = NULL,
                        ignore.reps.errors = FALSE) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # Libraries
  require(dplyr)
  require(ggplot2)

  # Reading the file
  if ("character" %in% class(results.file)) {
    results =
      readxl::read_excel(path = results.file,
                         sheet = results.sheet.position,
                         skip = rows.to.skip,
                         col_names = file.header)

    if (file.tail == TRUE) {# Removing last 5 rows
      results = results[-((nrow(results)-4):nrow(results)),]
    }
  } else {
    results = results.file
  }



  # Check reference sample
  if (!is.null(reference.sample) | (length(reference.sample) > 1)) {
    if (!(reference.sample %in% results$`Sample Name`)) {
      return(warning(paste0("The 'reference.sample' must be among the sample_IDs present in your file -> ",
                            paste(unique(results$`Sample Name`), collapse = ", "))))
    }
  } else {reference.sample = results$`Sample Name`[[1]]}



  # Check excluded samples
  if (!is.null(exclude.samples)) {
    if (F %in% (exclude.samples %in% results$`Sample Name`)) {
      message(paste0("The 'exclude.samples' values are not among the sample_IDs present in your file -> ",
                     paste(unique(results$`Sample Name`), collapse = ", "), ". \nHowever, the function will run anyway."))
    }
  }



  # Define the sample order
  actual_order = unique(results$`Sample Name`)

  if (is.null(samples.order)) {
    samples.order = c(reference.sample, actual_order[actual_order != reference.sample])
  } else {samples.order = samples.order}

  results =
    dplyr::mutate(.data = results, `Sample Name` = factor(`Sample Name`, levels = samples.order)) %>%
    dplyr::arrange(`Sample Name`, `Target Name`)



  # Set 'undetermined' samples to 0
  results =
    results %>%
    dplyr::mutate(CT = gsub(pattern = "undetermined" ,
                            replacement = 0,
                            x = results$CT,
                            ignore.case = T)) %>%
    dplyr::mutate(CT = as.numeric(CT))


  # Check housekeeping genes
  if (F %in% (housekeeping.genes %in% results$`Target Name`) | is.null(housekeeping.genes)) {
    return(warning(paste0("The 'housekeeping.genes' must be among the targets present in your file -> ",
                          paste(unique(results$`Target Name`), collapse = ", "))))
  }


  # Re-shape table for replicates
  reps_tb = data.frame()

  for (s in levels(results$`Sample Name`)) {
    for (t in unique(results$`Target Name`)) {
      current_tb = dplyr::filter(.data = results,
                                 `Sample Name` == s,
                                 `Target Name` == t)

      current_reps_tb = data.frame(t(data.frame(c(s, t, current_tb$CT))))
      colnames(current_reps_tb) = c("Sample Name", "Target Name", paste0("rep", 1:(ncol(current_reps_tb)-2)))

      reps_tb = plyr::rbind.fill(reps_tb, current_reps_tb)
    }
  }

  rownames(reps_tb) = NULL
  reps_tb[is.na(reps_tb)] = 0


  # Compute reps differences
  combinations = combn(colnames(reps_tb[-c(1:2)]), 2)

  for (c in 1:ncol(combinations)) {
    reps_tb =
      reps_tb %>%
      dplyr::mutate(diff = abs(as.numeric(reps_tb[,combinations[1,c]]) - as.numeric(reps_tb[,combinations[2,c]])))

    colnames(reps_tb)[ncol(reps_tb)] = paste0(combinations[1,c], "-", combinations[2,c])
  }


  # Remove too different CTs
  reps_tb_clean = reps_tb
  rep_names = unique(c(combinations[1, ], combinations[2, ]))
  to_keep_tb = reps_tb[(ncol(reps_tb) - ncol(combinations) + 1):ncol(reps_tb)] < max.delta.reps

  list_to_remove = list()
  for (i in 1:nrow(to_keep_tb)) {
    failed = sum(to_keep_tb[i, ] == F)

    if (failed != 0) {
      if ((failed == (ncol(combinations) - 2)) | (failed == ncol(combinations))) {
        list_to_remove[[i]] = c("all")
      }
      else {
        to_remove = names(to_keep_tb[i, to_keep_tb[i,] == FALSE])
        to_remove = strsplit(to_remove, "-")
        list_to_remove[[i]] = c(unique(unlist(to_remove)[duplicated(unlist(to_remove))]))
      }
    }
    else {
      if (failed == 0) {list_to_remove[[i]] = c("none")}
    }
  }



  # Substitute the reps-to-remove with NA
  for (i in 1:length(list_to_remove)) {
    if (list_to_remove[[i]][1] != "none") {
      if (list_to_remove[[i]][1] == "all") {
        for (k in 1:length(rep_names)) {
          reps_tb_clean[i, rep_names[k]] = NA
        }
      }
      else {
        for (k in 1:length(list_to_remove[[i]])) {
          reps_tb_clean[i, list_to_remove[[i]][k]] = NA
        }
      }
    }
  }



  # Plot the reps difference
  differences_table = data.frame()

  for (r in paste0(combinations[1,], "-", combinations[2,])) {
    differences_table = rbind(differences_table,
                              data.frame(Sample_Name = reps_tb$`Sample Name`,
                                         Target_Name = reps_tb$`Target Name`,
                                         delta = reps_tb[,r],
                                         comparison = r))
  }

  differences_table =
    differences_table %>%
    dplyr::mutate(Sample_Name = factor(Sample_Name, levels = rev(levels(results$`Sample Name`))),
                  good = factor(delta < max.delta.reps, levels = c(T, F)))


  reps_difference_plot =
    ggplot(differences_table,
           aes(y = Sample_Name,
               x = comparison)) +
    geom_tile(aes(fill = good),
              color = "white",
              size = 1.5) +
    geom_text(aes(label = round(delta, 2))) +
    scale_fill_manual(values = c("#f1f5ff", "#ff8e8e"), drop = F) +
    facet_wrap(~Target_Name) +
    theme_classic() +
    theme(axis.text = element_text(color = "#000000"),
          axis.text.x = element_text(angle = x.labels.rotation,
                                     hjust = 1))



  # Decide which table to use
  if (ignore.reps.errors == T) {analyses_tb = reps_tb} else {analyses_tb = reps_tb_clean}

  # Convert reps values to numeric
  for (i in rep_names) {analyses_tb[,i] = as.numeric(analyses_tb[,i])}

  # Compute means of replicates
  analyses_tb = dplyr::mutate(.data = analyses_tb,
                              CT_reps_mean = rowMeans(as.matrix(analyses_tb[,rep_names]), na.rm = T))




  ##########################
  # Compute the 2^-(target-housekeeping) and FoldChange
  list_norm_tb = list()

  for (h in housekeeping.genes) {
    current_normalized_analyses_tb = data.frame()

    for (t in unique(results$`Target Name`)) {
      # Subset the analyses table
      current_gene = analyses_tb %>% dplyr::filter(`Target Name` == t)

      current_house =
        analyses_tb %>%
        dplyr::filter(`Target Name` == h) %>%
        dplyr::mutate(house_mean = CT_reps_mean) %>%
        dplyr::select(`Sample Name`, house_mean)

      current_norm =
        dplyr::left_join(current_gene, current_house, by = "Sample Name") %>%
        dplyr::mutate(norm_exp = 2^-(CT_reps_mean-house_mean))

      current_norm =
        current_norm %>%
        dplyr::mutate(FoldChange = norm_exp / current_norm[which(current_norm$`Sample Name` == reference.sample),"norm_exp"]) %>%
        dplyr::select(-house_mean)

      colnames(current_norm)[(ncol(current_norm)-1):ncol(current_norm)] = c(paste0("exp_norm_to_",h),
                                                                            paste0("FC_over_", reference.sample))
      current_normalized_analyses_tb = rbind(current_normalized_analyses_tb, current_norm)
    }

    list_norm_tb[[which(housekeeping.genes == h)]] = current_normalized_analyses_tb
    names(list_norm_tb)[which(housekeeping.genes == h)] = h
  }




  #############################
  # Generate the plots
  exp_plots = list()
  FC_plots = list()

  for (h in housekeeping.genes) {
    tb =
      list_norm_tb[[h]] %>%
      dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = levels(results$`Sample Name`))) %>%
      dplyr::filter(!(`Sample Name` %in% exclude.samples))

    colnames(tb)[grepl("exp_norm_to_|FC_over_", colnames(tb))] = c("exp", "FC")


    exp_plots[[which(housekeeping.genes == h)]] =
      ggplot(data = tb %>% dplyr::filter(`Target Name` != h),
             aes(x = `Sample Name`,
                 y = exp)) +
      geom_bar(stat = "identity", fill = "steelblue", position = position_dodge(width=0.9)) +
      geom_text(aes(label=round(exp,3)), position = position_dodge(width=0.9), vjust=-0.25, size = text.size) +
      facet_wrap(~`Target Name`,
                 scales = ifelse(test = fix.y.axis == T,
                                 yes = "fixed",
                                 no = "free")) +
      ggtitle(paste0("Normalization relative to ", h)) +
      ylab("RNA expression") +
      theme_classic() +
      theme(axis.text = element_text(color = "#000000"),
            axis.text.x = element_text(angle = x.labels.rotation,
                                       hjust = 1))


    # Remove housekeeping.genes from FC plots (if required)
    if (exclude.housekeeping.FC == T) {tb = dplyr::filter(.data = tb, !(`Target Name` %in% housekeeping.genes))}


    FC_plots[[which(housekeeping.genes == h)]] =
      ggplot(data = tb,
             aes(x = `Sample Name`,
                 y = FC)) +
      geom_bar(stat = "identity", fill = "#ffb871", position = position_dodge(width=0.9)) +
      geom_text(aes(label=round(FC,2)), position = position_dodge(width=0.9), vjust=-0.25, size = text.size) +
      facet_wrap(~`Target Name`,
                 scales = ifelse(test = fix.y.axis == T,
                                 yes = "fixed",
                                 no = "free")) +
      ylab(paste0("FC relative to ", reference.sample)) +
      ggtitle(paste0("Normalization relative to ", h)) +
      theme_classic() +
      geom_hline(yintercept = 1, linetype = 2, color = "#000000") +
      theme(axis.text = element_text(color = "#000000"),
            axis.text.x = element_text(angle = x.labels.rotation,
                                       hjust = 1))


    names(exp_plots)[which(housekeeping.genes == h)] = h
    names(FC_plots)[which(housekeeping.genes == h)] = h
  }




  ###############
  # COMPUTE THE MEANS BETWEEN HOUSEKEEPING

  FC_matrix = list()

  for (h in 1:length(housekeeping.genes)) {
    FC_matrix[[h]] = list_norm_tb[[h]][which(grepl("FC_over", colnames(list_norm_tb[[h]])))]
  }

  FC_matrix =
    data.frame(FC_matrix) %>%
    dplyr::mutate(mean_FC_over_refereceSample = rowMeans(as.matrix(data.frame(FC_matrix)), na.rm = T),
                  mean_FC_SEM = matrixStats::rowSds(as.matrix(data.frame(FC_matrix)), na.rm = T) / sqrt(length(housekeeping.genes)))

  mean_table =
    list_norm_tb[[1]][,-((ncol(list_norm_tb[[1]])-3):ncol(list_norm_tb[[1]]))] %>%
    mutate(mean_FC_over_refereceSample = FC_matrix$ mean_FC_over_refereceSample,
           mean_FC_SEM = FC_matrix$mean_FC_SEM)

  names(mean_table)[which(colnames(mean_table) == "mean_FC_over_refereceSample")] = paste0("mean_FC_over_", reference.sample)

  list_norm_tb[[length(list_norm_tb)+1]] = mean_table
  names(list_norm_tb)[length(list_norm_tb)] = "mean_FC_housekeeping"


  # generate mean plots
  tb =
    list_norm_tb$mean_FC_housekeeping %>%
    dplyr::mutate(`Sample Name` = factor(`Sample Name`, levels = samples.order)) %>%
    dplyr::filter(!(`Sample Name` %in% exclude.samples),
                  !(`Target Name` %in% housekeeping.genes))

  colnames(tb)[grepl("FC_over_", colnames(tb))] = "FC"


  FC_plots[[length(FC_plots)+1]] =
    ggplot(data = tb,
           aes(x = `Sample Name`,
               y = FC)) +
    geom_bar(stat = "identity", fill = "#ff9225", position = position_dodge(width=0.9)) +
    geom_errorbar(aes(ymin = FC-mean_FC_SEM, ymax = FC+mean_FC_SEM),
                  position = position_dodge(width=0.9),
                  width = 0.3, color = "#000000") +
    geom_text(data = tb,
              aes(x = `Sample Name`,
                  y = FC+mean_FC_SEM,
                  label  =round(FC,2)),
              position = position_dodge(width=0.9),
              vjust=-0.25, size = text.size,
              inherit.aes = F) +
    facet_wrap(~`Target Name`,
               scales = ifelse(test = fix.y.axis == T,
                               yes = "fixed",
                               no = "free")) +
    ylab("mean FoldChange \U00B1 SEM") +
    ggtitle("Normalization relative to all housekeeping") +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = 2, color = "#000000") +
    theme(axis.text = element_text(color = "#000000"),
          axis.text.x = element_text(angle = x.labels.rotation,
                                     hjust = 1))

  names(FC_plots)[length(FC_plots)] = "mean.FC.housekeeping"



  #########################
  # Export the output

  return(list(original.table = results,
              reshaped.table = reps_tb,
              reshaped.table.cleaned = reps_tb_clean,
              reps.validation.plot = reps_difference_plot,
              analyzed.data = list_norm_tb,
              expression.plots = exp_plots,
              foldChange.plots = FC_plots))

} # END function
