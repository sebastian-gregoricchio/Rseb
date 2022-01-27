#' @title  Single base bigWig score selector
#'
#' @description Function to get the score from a bigWig for each base in a given genomic region.
#'
#' @param region An atomic string indicating the genomic region into which restrict the final plot in the format 'chr1:1234-5678'.
#' @param bigWig Full path to a bigWig file.
#' @param missing.data.as.zero A logical value to define whether missing data (NAs) should be treated as zeros. By default \code{TRUE}.
#' @param reverse.score A logical value to indicate whether the score order should be inverted. Default \code{TRUE}.
#'
#' @return The output is a numeric vector containing the score for each base at a given position.
#'
#' @export get.single.base.score.bw
#'
#'
#'
get.single.base.score.bw =
  function(region,
           bigWig,
           missing.data.as.zero = TRUE,
           reverse.score = FALSE) {

    # -------------------------------------------------------------------------------------
    ### LIBRARIES
    require(rtracklayer)
    require(GenomicRanges)
    require(IRanges)
    # -------------------------------------------------------------------------------------


    ### Check genomic region and generate a bed file
    if (class(region) != "character") {return(warning("The genomic region must be a string in the format 'chr1:12345-67890', where the end position must be greater than the start position and the window must have a length > 1."))}
    region = gsub(pattern = "[,]", replacement = "", x = region)
    region = gsub(pattern = "[-]|[+]|[:]", replacement = "_", x = region)
    region = strsplit(region, "_")[[1]]
    if (region[3] <= region[2]) {return(warning("The genomic region must be a string in the format 'chr1:12345-67890', where the end position must be greater than the start position and the window must have a length > 1."))}
    region = GRanges(seqnames = region[1], ranges = IRanges(start = as.numeric(region[2]), end = as.numeric(region[3])))



    ### get score
    # Check if the bigWig chrom names are the same than in the bed regions
    if (suppressWarnings(inherits(try(import(BigWigFile(bigWig), selection = region, as = 'NumericList')[[1]],
                                      silent = TRUE),
                                  "try-error"))) { # it is TRUE when it does not work due to chromosome names in the bigWig != bed
      region = diffloop::rmchr(region)
    }

    score = import(BigWigFile(bigWig), selection = region, as = 'NumericList')[[1]]
    if (reverse.score == TRUE) {score = rev(score)}

    if (missing.data.as.zero == TRUE) {score[is.na(score)]=0}

    return(score)

} # END FUNCTION
