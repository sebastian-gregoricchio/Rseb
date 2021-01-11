#' @title Conversion of GSEA terms into Gene Ontology numbers
#'
#' @description Helps to convert the terms of GSEA analyses into Gene Ontology (GO) ID numbers.
#'
#' @param input_terms A character vector containing the GSEA terms to be converted.
#' @param input_pvalue A numeric vector containing the p-values of the GSEA terms.
#' @param return_table Logic value to define whether to return the resulting data.frame. By default \code{TRUE}.
#' @param export_table Logic value to define whether to export the resulting data.frame. By default \code{FALSE}.
#' @param output_file_name Path and file name of the output table if export is required. By default \code{<working.directory>/GO_numbers_table.tsv}.
#'
#' @return If required, returns a data.frame with 3 columns: GO_number, GO_annotation, p.value. This table could be directly exported.
#'
#' @details This functions requires the package \code{GO.db}. \cr If problems are encountered during the installation see \url{https://www.biostars.org/p/50564/}.
#'
#' @export GSEA.to.GOnumber
#'
# @import GO.db
# @import dplyr


GSEA.to.GOnumber = function(input_terms,
                            input_pvalue,
                            return_table = T,
                            export_table = F,
                            output_file_name = paste(getwd(), "GO_numbers_table.tsv", sep ="/")) {

# Let's check if you have the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

  pkg = "GO.db"
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    if(!require(pkg, character.only = TRUE)) stop(paste(pkg,"package not found."))}

if (!require("dplyr",character.only = TRUE)) {
  install.packages("dplyr", dependencies = TRUE)
  if(!require("dplyr",character.only = TRUE)) stop("Package not found")}

#Retrieve the DB
go_DB = as.data.frame(Term(GOTERM))
go_DB = data.frame(GO_number = rownames(go_DB), GO_annotation = go_DB$`Term(GOTERM)`)

input_table = data.frame(terms = input_terms, pvalue = input_pvalue)

go_DB$GO_annotation = tolower(gsub(x = go_DB$GO_annotation, pattern = "-", replacement = " "))
go_DB$GO_annotation = tolower(gsub(x = go_DB$GO_annotation, pattern = ",", replacement = ""))
input_table$terms = tolower(gsub(x = input_table$terms, pattern = "GO_", replacement = ""))
input_table$terms = tolower(gsub(x = input_table$terms, pattern = "_", replacement = " "))

go_numbers = go_DB %>% filter(go_DB$GO_annotation %in% input_table$terms)
go_numbers = go_numbers %>% arrange(go_numbers$GO_annotation)
pvalues_column = input_table %>% filter(input_table$terms %in% go_numbers$GO_annotation)
pvalues_column = pvalues_column %>% arrange(pvalues_column$terms)
pvalues_column = pvalues_column$pvalue
go_numbers = go_numbers %>% mutate(p.value = pvalues_column)

if (export_table == T) {
  write.table(x = go_numbers,
              file = output_file_name,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE, col.names = TRUE)
}

if (return_table == T) {return(go_numbers)}

}
