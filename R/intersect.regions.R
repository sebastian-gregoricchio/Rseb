#' @title Genomic regions overlapper
#'
#' @description A tool to define overlaps between bed files/regions derived from different formats. The function allows the overlap in stranded mode and can considered a specific minimal percentage of overlap between regions.
#'
#' @param reference.regions A single value or a list of regions to be used as 'reference'. The values accepted are: a. a character with the full path to a bed file, b. a data.frame in at least BED3 format, c. a GRanges object in at least BED3 format. If a list of elements is provided all the regions will be merged in a unique combined list and only completely identical regions will be remove to avoid duplicates. Combination of different formats is allowed.
#' @param test.regions A single value or a list of regions to be used as 'test'. The values accepted are: a. a character with the full path to a bed file, b. a data.frame in at least BED3 format, c. a GRanges object in at least BED3 format. If a list of elements is provided all the regions will be merged in a unique combined list and only completely identical regions will be remove to avoid duplicates. Combination of different formats is allowed.
#' @param min.percentage.reference A numeric value in 0-100 to define which percentage of a region in the 'reference' dataset must overlap with a region in the 'test' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.percentage.test A numeric value in 0-100 to define which percentage of a region in the 'test' dataset must overlap with a region in the 'reference' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.bases.overlap Integer, greater than 0, value to indicate the minimal number of bases to consider as minimum overlap between two regions. Non integer values will be rounded at integer, while number lower that 1 will be coerced to 1. Default value: \code{1}.
#' @param sort.overlaps Logic value to define whether the output should be sorted or not. Default value: \code{FALSE}.
#' @param stranded A logical value to define whether the analyses should be performed by strand: regions in one strand will be overlapped only with regions of the same strand. The strand symbols considered are '+' and '-', any other symbol will considered in a unique separated category. Default value: \code{FALSE}.
#' @param return.as.data.frame Logical value to define whether the output list should contain data.frames instead of GRanges objects. Default value: \code{TRUE}.
#'
#' @return The function returns a list of data.frames/GRanges objects containing:
#' \itemize{
#'   \item \code{overlaps.reference}: XX;
#'   \item \code{non.overlaps.reference}: XXX;
#'   \item \code{overlaps.testt}: VV;
#'   \item \code{non.overlaps.test}: XX.
#'  }
#'
#' @export intersect.regions

# @import GenomicRanges
# @import diffloop
# @import dplyr
# @import IRanges
# @import S4Vectors

intersect.regions =
  function(reference.regions,
           test.regions,
           min.percentage.reference = 0,
           min.percentage.test = 0,
           min.bases.overlap = 1,
           sort.overlaps = FALSE,
           stranded = FALSE,
           return.as.data.frame = TRUE) {

    ######################################################################################
    ### Required libraries
    # require("GenomicRanges")
    # require(diffloop)
    require(dplyr)
    require(S4Vectors)
    require(GenomicRanges)
    require(IRanges)

    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)
    ######################################################################################


    ### Check of the thresholds
    if (min.percentage.reference < 0) {min.percentage.reference = 0}
    if (min.percentage.reference > 100) {min.percentage.reference = 100}

    if (min.percentage.test < 0) {min.percentage.test = 0}
    if (min.percentage.test > 100) {min.percentage.test = 100}

    min.bases.overlap = round(min.bases.overlap, 0)
    if (min.bases.overlap < 1) {min.bases.overlap = 1}


    ### loading regions and convert to GRanges object
    # Defining import/reading function
    read.regions = function(x) {
      if ("data.frame" %in% class(x)) {
        x = GenomicRanges::makeGRangesFromDataFrame(df = as.data.frame(x),
                                                    seqnames.field = colnames(x)[1],
                                                    start.field = colnames(x)[2],
                                                    end.field = colnames(x)[3],
                                                    strand.field = ifelse(test = ncol(x)>=6,
                                                                          yes = colnames(x)[6],
                                                                          no = "*"))
      } else if ("character" %in% class(x)) {
        x = rtracklayer::import.bed(x)
      } else if ("GRanges" %in% class(x)) {
        x = x
      } else {
        return(return(warning("The format of One of the regions provided is a non recognized. Formats allowed:\n   - 'data.frame' in at least BED3 format;\n   - 'characther' string with the full path to a .bed file;\n   - 'GRanges' bed object.")))
      }
      return(unique(diffloop::addchr(diffloop::rmchr(x))))
    } # end read.regions


    # Reading reference regions
    if ("list" %in% class(reference.regions)) {
      reference.regions = lapply(reference.regions, function(x){return(read.regions(x))})
      reference.regions = unique(Reduce(x = reference.regions, append))
    } else {
      reference.regions = read.regions(reference.regions)
    }


    # Reading test regions
    if ("list" %in% class(test.regions)) {
      test.regions = lapply(test.regions, function(x){return(read.regions(x))})
      test.regions = unique(Reduce(x = test.regions, append))
    } else {
      test.regions = read.regions(test.regions)
    }



    ### Performing the overlaps
    # Defining the function to perform the overlaps
    find.intersections =
      function(reference.regions,
               test.regions,
               min.percentage.reference,
               min.percentage.test,
               min.bases.overlap,
               sort.overlaps) {

        # Defining the overlaps in each data.set
        total.hits = IRanges::findOverlaps(query = reference.regions,
                                           subject = test.regions,
                                           minoverlap = min.bases.overlap)

        # Get only the overhanging regions
        overlaps = IRanges::pintersect(reference.regions[queryHits(total.hits)],
                                       test.regions[subjectHits(total.hits)])

        # Compute the fraction of overlap for a
        overalp.perc =
          cbind(data.frame(overlaps),
                data.frame(total.hits) %>%
                  dplyr::mutate(perc.overlap.reference = (GenomicRanges::width(overlaps) / GenomicRanges::width(reference.regions[queryHits(total.hits)]))*100,
                                perc.overlap.test = (GenomicRanges::width(overlaps) / GenomicRanges::width(test.regions[subjectHits(total.hits)]))*100))



        filtered.hits = total.hits[overalp.perc$perc.overlap.reference >= min.percentage.reference &
                                     overalp.perc$perc.overlap.test >= min.percentage.test]



        # Return Overlaps in ref, overlaps in test, and the regions not overlapping for each region list
        overlaps.reference = unique(reference.regions[queryHits(filtered.hits)])
        overlaps.test = unique(test.regions[subjectHits(filtered.hits)])

        overlaps.list = list(overlaps.reference = overlaps.reference,
                             non.overlaps.reference = unique(reference.regions[!(reference.regions %in% overlaps.reference)]),
                             overlaps.test = overlaps.test,
                             non.overlaps.test = unique(test.regions[!(test.regions %in% overlaps.test)]))

        if (isTRUE(sort.overlaps)) {
          overlaps.list$overlaps.reference = sort(overlaps.list$overlaps.reference)
          overlaps.list$non.overlaps.reference = sort(overlaps.list$non.overlaps.reference)
          overlaps.list$overlaps.test = sort(overlaps.list$overlaps.test)
          overlaps.list$non.overlaps.test = sort(overlaps.list$non.overlaps.test)
        }

        return(overlaps.list)
      } # end find.intersections



    # Find the intersections
    if (isTRUE(stranded)) {
      # Check that only at least the "+" strand is present
      if (!(("+" %in% reference.regions@strand | "-" %in% reference.regions@strand) & ("+" %in% test.regions@strand | "-" %in% test.regions@strand))) {
        message("None of the regions provided contain any '+' or '-' symbol in the strand column, therefore the 'stranded' option will be ignored.")
      } else {
        # split Granges by strand
        reference.plus = reference.regions[reference.regions@strand == "+"]
        reference.minus = reference.regions[reference.regions@strand == "-"]
        reference.other = reference.regions[!(reference.regions@strand %in% c("+","-"))]

        test.plus = test.regions[test.regions@strand == "+"]
        test.minus = test.regions[test.regions@strand == "-"]
        test.other = test.regions[!(test.regions@strand %in% c("+","-"))]


        # find overlaps by strand
        overlaps.plus = find.intersections(reference.regions = reference.plus,
                                           test.regions = test.plus,
                                           min.percentage.reference = min.percentage.reference,
                                           min.percentage.test = min.percentage.test,
                                           min.bases.overlap = min.bases.overlap,
                                           sort.overlaps = sort.overlaps)

        overlaps.minus = find.intersections(reference.regions = reference.minus,
                                            test.regions = test.minus,
                                            min.percentage.reference = min.percentage.reference,
                                            min.percentage.test = min.percentage.test,
                                            min.bases.overlap = min.bases.overlap,
                                            sort.overlaps = sort.overlaps)

        overlaps.other = find.intersections(reference.regions = reference.other,
                                            test.regions = test.other,
                                            min.percentage.reference = min.percentage.reference,
                                            min.percentage.test = min.percentage.test,
                                            min.bases.overlap = min.bases.overlap,
                                            sort.overlaps = sort.overlaps)

        final.overlaps = list(overlaps.reference = unique(Reduce(x = list(overlaps.plus$overlaps.reference,
                                                                          overlaps.minus$overlaps.reference,
                                                                          overlaps.other$overlaps.reference),
                                                                 append)),
                              non.overlaps.reference = unique(Reduce(x = list(overlaps.plus$non.overlaps.reference,
                                                                              overlaps.minus$non.overlaps.reference,
                                                                              overlaps.other$non.overlaps.reference),
                                                                     append)),
                              overlaps.test = unique(Reduce(x = list(overlaps.plus$overlaps.test,
                                                                     overlaps.minus$overlaps.test,
                                                                     overlaps.other$overlaps.test),
                                                            append)),
                              non.overlaps.test = unique(Reduce(x = list(overlaps.plus$non.overlaps.test,
                                                                         overlaps.minus$non.overlaps.test,
                                                                         overlaps.other$non.overlaps.test),
                                                                append)))
      } # end stranded
    } else {
      final.overlaps = find.intersections(reference.regions = reference.regions,
                                          test.regions = test.regions,
                                          min.percentage.reference = min.percentage.reference,
                                          min.percentage.test = min.percentage.test,
                                          min.bases.overlap = min.bases.overlap,
                                          sort.overlaps = sort.overlaps)
    }



    ### Return the overlaps
    if (isTRUE(return.as.data.frame)) {
      return(list(overlaps.reference = as.data.frame(final.overlaps$overlaps.reference),
                  non.overlaps.reference = as.data.frame(final.overlaps$non.overlaps.reference),
                  overlaps.test = as.data.frame(final.overlaps$overlaps.test),
                  non.overlaps.test = as.data.frame(final.overlaps$non.overlaps.test)))
    } else {
      return(final.overlaps)
    }
  } # END of the function
