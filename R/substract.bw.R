#' @title Combination of two or more list in a unique one.
#'
#' @description Combines two or more lists in a single one keeping the elements names
#'
#' @param bw1 Full path to the first bigWig (the second one will be substracted to this one).
#' @param bw2 Full path to the second bigWig (it will be substracted to the first one).
#' @param return.substracted.bw Logic value to define whether return the resulting bigWig as GRanges object. By default \code{FALSE}.
#' @param export.substracted.bw Logic value to define whether export the resulting bigWig. By default \code{TRUE}.
#' @param substracted.bw.file String for the path of the resulting bigwig file to be exported. \cr By default \code{<working.directory>/subtraction.bw}.
#'
#' @return If required a subtraction bigWig is returned as GRanges object. The resulting bigWig can be also directly exported.
#'
#' @export substract.bw
#'
#' @import rtracklayer


substract.bw = function(bw1,
                        bw2,
                        wd = getwd(),
                        return.substracted.bw = F,
                        export.substracted.bw = T,
                        substracted.bw.file = paste(getwd(), "subtraction.bw", sep = "/")) {

  # Loading libraries
  #require(sf) # due to a bug load this before
  pkg = "rtracklayer"
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    if(!require(pkg, character.only = TRUE)) stop(paste(pkg,"package not found."))}

  require("rtracklayer")

  # Import/read the BigWigs
  bws = list(bw1, bw2)

  for (i in 1:length(bws)) {
    if (class(bws[[i]]) != "GRanges" & class(bws[[i]]) == "character") {
      bws[[i]] = import(bws[[i]])} else if (class(bws[[i]]) != "GRanges" & class(bws[[i]]) != "character") {
        return(warning("Invalid input format. Either GRanges bigWig or a string with the path for the file"))
      }
  }


  # Taking only the overlapping regions
  if (length(bw1) != length(bw2)) {
    if (length(bw1) > length(bw2)) {
      message("Overlapped subgrouping ...")
      bw1 = bw1[bw1 %in% bw2]
    } else {
      message("Overlapped subgrouping ...")
      bw2 = bw2[bw2 %in% bw1]
    }
  }


  # Subtraction of scores >>> bw1 - bw2
  score(bws[[1]]) = score(bws[[1]]) - score(bws[[2]])

  # Export the subtracted file if required
  if (export.substracted.bw == T) {
    export(bws[[1]], substracted.bw.file)
    message("The result was exported has", substracted.bw.file)
  }

  # Return the subtracted bigWig if required
  if (return.substracted.bw == T) {
    return(bws[[1]])
  }

} # end function
