#' @title Conversion of ENSEMBL gene IDs.
#'
#' @description Conversion of ENSEMBL gene IDs to gene symbols.
#'
#' @param ensembl.id String vector of ENSEMBL genes IDs
#' @param type String to define the type of ENSEMBL inputs. By default \code{gene} to indicate "ensembl_gene_id". If different from "gene" it will be set to "ensembl_transcript_id_version".
#' @param organism String to define de organism, e.g. \code{mmusculus}, \code{hsapiens}, etc. By default \code{mmusculus}.
#'
#' @return A string vector with the corresponding gene_symbols.
#'
#' @examples
#' gene_symbols =
#' get.gene.name(
#'    ensembl.id = c("ENSMUSG00000002111", "ENSMUSG00000027381"),
#'    type = "gene",
#'    organism = "mmusculus")
#'
#' @export get.gene.name

get.gene.name = function(ensembl.id,
                         type = "gene",
                         organism = "mmusculus") { # gene or transcript

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  require(biomaRt)

  dataset = paste(organism, "_gene_ensembl", sep = "")

  mart = useMart(dataset = dataset,
                 biomart = "ensembl")

  filter = ifelse(test = type == "gene",
                  yes = "ensembl_gene_id",
                  no = "ensembl_transcript_id_version")

  gene_symbol = getBM(attributes = "external_gene_name",
                      filters = filter,
                      values = ensembl.id,
                      mart = mart)

  return(gene_symbol)

  } # end function
