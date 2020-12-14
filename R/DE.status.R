#' @title Differential Expression status calculator for RNA-seq data
#'
#' @description Defines the differential expression status of genes from RNA-seq data depending on fold change expression and adjusted p-value.
#'
#' @param log2FC Numeric vector of log2(fold change expression) values.
#' @param p.value.adjusted Numeric vector of p-values. Use of adjusted p-values is recommended.
#' @param FC_threshold Value of the threshold to use for the fold change expression to define differentially expressed genes, expressed as linear value. By default 1.5 and by consequence 1/1.5.
#' @param FC_NoResp_left Value of the threshold to use for the fold change expression to define unresponsive genes when \code{FC < 1}, expressed as linear value. By default 0.9. If \code{NULL} it will be calculated symmetrically from \code{FC_NoResp_rigth} as 1/\code{FC_NoResp_rigth}.
#' @param FC_NoResp_rigth Value of the threshold to use for the fold change expression to define unresponsive genes when \code{FC > 1}, expressed as linear value. By default 1.1. If \code{NULL} it will be calculated symmetrically from \code{FC_NoResp_left} as 1/\code{FC_NoResp_left}.
#' @param p.value_threshold Value of the threshold to use for the p-values to define differentially expressed genes, expressed as linear value. By default 0.05.
#' @param low.FC.status.label String to define the label indicating the differentially expressed genes with a \code{FoldChange < FC_threshold}.
#' @param high.FC.status.label String to define the label indicating the differentially expressed genes with a \code{FoldChange > FC_threshold}.
#' @param unresponsive.label String to define the label indicating the unresponsive genes identified as \code{FC_NoResp_left < FoldChange < FC_NoResp_rigth} and \code{p.value > p.value.threshold}.
#' @param null.label String to define the label indicating the null genes.
#'
#' @return It returns a vector containing the differential expression status for each original value in the same order used in the input.
#'
#' @export DE.status


DE.status = function(log2FC, # log2(FC)
                     p.value.adjusted,
                     FC_threshold = 1.5,
                     FC_NoResp_left = 0.9,
                     FC_NoResp_rigth = NULL,
                     p.value_threshold = 0.05,
                     low.FC.status.label = "DOWN",
                     high.FC.status.label = "UP",
                     unresponsive.label = "NoResp",
                     null.label = "NULL") {

  # Check the unresponsive thresholds
  if (is.null(FC_NoResp_left) & is.null(FC_NoResp_rigth)) {
    return(warning("'FC_NoResp_left' and 'FC_NoResp_rigth' parameters can't be both NULL."))
  }


  # Calculating the thresholds in log
  Fc = log2(FC_threshold)
  FcNS_left = log2(FC_NoResp_left)

  if (is.null(FC_NoResp_rigth)) {
    FcNS_rigth = -FcNS_left} else {
      FcNS_rigth = log2(FC_NoResp_rigth)}

  # Build a data.frame containing the input and an empty column for the DE status
  df = data.frame(FC_col = log2FC,
                  padj_col = p.value.adjusted,
                  UP.DOWN = null.label,
                  stringsAsFactors = F)

  # Definition of the 4 status: UP, DOWN, NoResp, NULL
  for (i in 1:nrow(df)) {
    if (df$FC_col[i] >= Fc & df$padj_col[i] <= p.value_threshold & !(is.na(df$padj_col[i]))) {
      df$UP.DOWN[i] = high.FC.status.label} else {

        if (df$FC_col[i] <= -Fc & df$padj_col[i] <= p.value_threshold & !(is.na(df$padj_col[i]))) {
          df$UP.DOWN[i] = low.FC.status.label} else {

            if ((FcNS_left <= df$FC_col[i]) & (df$FC_col[i] <= FcNS_rigth) & ((df$padj_col[i] > p.value_threshold) | (is.na(df$padj_col[i])))) {
              df$UP.DOWN[i] = unresponsive.label} else {

                df$UP.DOWN[i] = null.label}
          }
      }
  } # end for loop

  # Returns the vector of the DE status
  return(df$UP.DOWN)}
