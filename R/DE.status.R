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

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # Check the unresponsive thresholds
  if (is.null(FC_NoResp_left) & is.null(FC_NoResp_rigth)) {
    return(warning("'FC_NoResp_left' and 'FC_NoResp_rigth' parameters can't be both NULL."))
  }


  # Calculating the thresholds in log
  FC.th = log2(FC_threshold)
  FcNS_left = log2(FC_NoResp_left)

  if (is.null(FC_NoResp_rigth)) {
    FcNS_rigth = -FcNS_left} else {
      FcNS_rigth = log2(FC_NoResp_rigth)}


  # Define signif status function
  status =
    function(FC, p) {
      ifelse(p < p.value_threshold,
             yes = ifelse(abs(FC) >= FC.th,
                          yes = ifelse(sign(FC) == 1,
                                       yes = high.FC.status.label,
                                       no = low.FC.status.label),
                          no = null.label),
             no = ifelse(FC >= FcNS_left & FC <= FcNS_rigth,
                         yes = unresponsive.label,
                         no = null.label))
    }

  # Define the label for each FC x Padj combination
  diff.status =
    unlist(purrr::pmap(.l = list(FC = log2FC,
                                 p = p.value.adjusted),
                       .f = function(FC,p){status(FC,p)}))

  # Returns the vector of the DE status
  return(diff.status)
} # END function
