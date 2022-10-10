#' @title VennDiagram from region overlaps
#'
#' @description A tool to plot VennDiagrams from overlaps between bed files/regions derived from different formats. The function allows the overlap in stranded mode and can considered a specific minimal percentage of overlap between regions.
#'
#' @param region.list A list of regions to be used as to compute the overlaps. The values accepted are: a. a character with the full path to a bed file, b. a data.frame in at least BED3 format, c. a GRanges object in at least BED3 format. If a list of elements is provided all the regions will be merged in a unique combined list and only completely identical regions will be remove to avoid duplicates. Combination of different formats is allowed.
#' @param region.names String vector with the names of the regions in the order.
#' @param colors Vector to define the line and error area colors. If only one value is provided it will applied to all the samples/groups. If the number of values is lower than the the required one, a random set of colors will be generated. All standard R.colors values are accepted. By default \code{c("#00A5CF", "#F8766D", "#AC88FF", "#E08B00", "#00BA38", "#BB9D00", "#FF61C9", "gray30")}.
#' @param color.transparency Numeric floating value between 0-1 to indicate the transparency, aka alpha, of the colors.
#' @param min.percentage.reference A numeric value in 0-100 to define which percentage of a region in the 'reference' dataset must overlap with a region in the 'test' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.percentage.test Numeric value in 0-100 to define which percentage of a region in the 'test' dataset must overlap with a region in the 'reference' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.
#' @param min.bases.overlap Integer, greater than 0, value to indicate the minimal number of bases to consider as minimum overlap between two regions. Non integer values will be rounded at integer, while number lower that 1 will be coerced to 1. Default value: \code{1}.
#' @param input.type String with the type of input provided to the euler function. Available values are \code{union} and \code{disjoint}. Default value: \code{union}.
#' @param shape.type String with the type of shape to use for the plot: one among \code{ellipse} and \code{circle}. Default value: \code{ellipse}.
#' @param plot.quantities Logical value to indicate whether the quantity of each subintersection should be plotted or not. By default \code{TRUE}.
#' @param stranded Logical value to define whether the analyses should be performed by strand: regions in one strand will be overlapped only with regions of the same strand. The strand symbols considered are '+' and '-', any other symbol will considered in a unique separated category. Default value: \code{FALSE}.
#'
#' @return The output is the Venn Diagram in an object of class eulergram/gTree/grob/gDesc.
#'
#' @export venn.overlap

# @import eulerr
# @import diffloop
# @import pryr


venn.overlap =
  function(region.list,
           region.names = LETTERS[1:length(region.list)],
           colors = c("#00A5CF", "#F8766D", "#AC88FF", "#E08B00", "#00BA38", "#BB9D00", "#FF61C9", "gray30"),
           color.transparency = 0.25,
           min.percentage.reference = 0,
           min.precentage.test = 0,
           min.bases.overlap = 1,
           input.type = "union",
           shape.type = "ellipse",
           plot.quantities = TRUE,
           stranded = FALSE) {


    ######################################################################################
    ### Required libraries
    require(pryr)

    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)
    ######################################################################################



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
    if ("list" %in% class(region.list)) {
      region.list = lapply(region.list, function(x){return(read.regions(x))})
    } else {
      return(warning("The 'region.list' must be a list."))
    }



    # Assign names to the region
    if (length(region.list) <= length(region.names)) {
      names(region.list) = region.names
    } else {
      region.names = LETTERS[1:length(region.list)]
      message("The region.names are less than the regions provided. Names are coerced to capital letters in alphabetic order.")
    }


    # Select colors
    if (length(colors) > length(region.list)) {
      colors = c("#00A5CF", "#F8766D", "#AC88FF", "#E08B00", "#00BA38", "#BB9D00", "#FF61C9", "gray30")[1:length(region.list)]
    } else {
      if (length(colors) < length(region.list)) {
        warning("The number of colors is not sufficient, random colors will be used.")

        colors = c()
        for (i in 1:length(region.list)){
          set.seed(runif(1,1,100))
          colors[i] = rgb(red = runif(1,0,1), green = runif(1,0,1), blue = runif(1,0,1))
        }
      }
    }




    # Compute all the combinations
    comb.list = list()
    for (i in 1:length(region.names)) {comb.list[[i]] = combn(x = region.names, m = i)}


    # Get the number of overlaps for each combination
    overalps = c()
    idx = 0

    for (c in 1:length(comb.list)) {
      comb.group = comb.list[[c]]

      for (k in 1:ncol(comb.group)) {
        idx = idx+1
        if (length(comb.group[,k]) == 1) {
          overalps[idx] = nrow(data.frame(region.list[[comb.group[,k]]]))
          names(overalps)[idx] = paste0(comb.group[,k], sep="", collapse = "&")
        } else {
          current.overlap = region.list[[comb.group[1,k]]]
          for (o in 2:length(comb.group[,k])) {
            current.overlap = (Rseb::intersect.regions(reference.regions = current.overlap,
                                                       test.regions = region.list[[comb.group[o,k]]],
                                                       min.percentage.reference = min.percentage.reference,
                                                       min.percentage.test = min.precentage.test,
                                                       min.bases.overlap = min.bases.overlap,
                                                       stranded = stranded,
                                                       return.as.data.frame = F,
                                                       sort.overlaps = F))$overlaps.reference
          }
          overalps[idx] = nrow(data.frame(current.overlap))
          names(overalps)[idx] = paste0(comb.group[,k], sep="", collapse = "&")
        }
      }
    }


    # Generate the distance/size matrix
    euler.mat = eulerr::euler(overalps, input = input.type, shape.type = shape.type)

    # Plotting the venn
    venn.plot %<a-% plot(euler.mat,
                         quantities = plot.quantities,
                         fill = colors,
                         alpha = color.transparency)


    # Return the plot
    return(venn.plot)

  } # END function
