#' @title Sorter function for .bed files.
#'
#' @description Sorts .bed files by chromosome and position.
#'
#' @param bed Two options are possible: \cr - String with the path to a .bed file; \cr - Data.frame corresponding to a bed file format (all the columns and their names will be kept).
#' @param bed.header Logic value to define whether the .bed file contains an header or not. By default \code{FALSE}.
#' @param sep String containing the separator character for a .bed file. By default \code{"\t"}.
#' @param return.bed Logic value to define if to return the bed as a data.frame. By default \code{TRUE}.  Only unique rows are kept.
#' @param export.bed Logic value to define if to export the bed file. By default \code{FALSE}. Only unique rows are kept.
#' @param export.file.name String to define the path to the file to be exported, if required. By default \code{"<working.directory>/sorted.bed"}.
#' @param export.header Logic value to define whether the header should be exported in the sorted bed file. By default \code{FALSE}.
#'
#' @return If required, returns a data.frame corresponding to the sorted .bed file.
#'
#' @details The function keeps only unique rows. \cr To get more information about the bed file format see the following page: \cr \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}.
#'
#' @export sort.bed
#'
# @import dplyr


sort.bed = function(bed,
                    bed.header = F,
                    sep = "\t",
                    return.bed = T,
                    export.bed = F,
                    export.file.name = paste(getwd(), "sorted.bed", sep="/"),
                    export.header = F) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  require(dplyr)

  if (class(bed)[1] == "character") {
    bed = read.delim(file = bed,
                     header = bed.header,
                     sep = sep,
                     stringsAsFactors = F)
    }

  # Sort bed for chromosomes 1-9
  bed.1_9 =
    bed %>%
    filter(grepl(pattern = "chr[0-9]$",
                 x = bed[,1]))
  bed.1_9 =
    bed.1_9 %>%
    arrange(bed.1_9[,1], bed.1_9[,2], bed.1_9[,3])


  # Sort bed for chromosomes 10-10+
  bed.10_more =
    bed %>%
    filter(grepl(pattern = "chr[0-9][0-9]$",
                x = bed[,1]))
  bed.10_more =
    bed.10_more %>%
    arrange(bed.10_more[,1], bed.10_more[,2], bed.10_more[,3])


  # Sort bed for non numeric chromosomes
  bed.others =
    bed %>%
    filter(!grepl(pattern = "chr[0-9]$|chr[0-9][0-9]$",
                  x = bed[,1]))
  bed.others =
    bed.others %>%
    arrange(bed.others[,1], bed.others[,2], bed.others[,3])


  # Merge the all sub beds
  sorted.bed = rbind(bed.1_9, bed.10_more, bed.others)


  # Check rows
  if (nrow(sorted.bed) != nrow(bed)) {
    return(warning("Something went wrong: the number of rows of the original bed is not the same of the sorted file."))
  }

  # Export if required
  if (export.bed == T) {
    write.table(x = unique(sorted.bed),
                file = export.file.name,
                quote = F, sep = sep,
                row.names = F, col.names = export.header)

    message(paste("The sorted file have been exported to",
                  export.file.name))
  }

  message("Only unique values are kept.")

  # Return the bed
  if (return.bed == T) {return(unique(sorted.bed))}
}
