% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venn.overlap.R
\name{venn.overlap}
\alias{venn.overlap}
\title{VennDiagram from region overlaps}
\usage{
venn.overlap(
  region.list,
  region.names = LETTERS[1:length(region.list)],
  colors = c("#00A5CF", "#F8766D", "#AC88FF", "#E08B00", "#00BA38", "#BB9D00", "#FF61C9",
    "gray30"),
  color.transparency = 0.25,
  min.percentage.reference = 0,
  min.percentage.test = 0,
  min.bases.overlap = 1,
  input.type = "union",
  shape.type = "ellipse",
  plot.quantities = TRUE,
  stranded = FALSE
)
}
\arguments{
\item{region.list}{A list of regions to be used as to compute the overlaps. The values accepted are: a. a character with the full path to a bed file, b. a data.frame in at least BED3 format, c. a GRanges object in at least BED3 format. If a list of elements is provided all the regions will be merged in a unique combined list and only completely identical regions will be remove to avoid duplicates. Combination of different formats is allowed.}

\item{region.names}{String vector with the names of the regions in the order.}

\item{colors}{Vector to define the line and error area colors. If only one value is provided it will applied to all the samples/groups. If the number of values is lower than the the required one, a random set of colors will be generated. All standard R.colors values are accepted. By default \code{c("#00A5CF", "#F8766D", "#AC88FF", "#E08B00", "#00BA38", "#BB9D00", "#FF61C9", "gray30")}.}

\item{color.transparency}{Numeric floating value between 0-1 to indicate the transparency, aka alpha, of the colors.}

\item{min.percentage.reference}{A numeric value in 0-100 to define which percentage of a region in the 'reference' dataset must overlap with a region in the 'test' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.}

\item{min.percentage.test}{Numeric value in 0-100 to define which percentage of a region in the 'test' dataset must overlap with a region in the 'reference' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.}

\item{min.bases.overlap}{Integer, greater than 0, value to indicate the minimal number of bases to consider as minimum overlap between two regions. Non integer values will be rounded at integer, while number lower that 1 will be coerced to 1. Default value: \code{1}.}

\item{input.type}{String with the type of input provided to the euler function. Available values are \code{union} and \code{disjoint}. Default value: \code{union}.}

\item{shape.type}{String with the type of shape to use for the plot: one among \code{ellipse} and \code{circle}. Default value: \code{ellipse}.}

\item{plot.quantities}{Logical value to indicate whether the quantity of each subintersection should be plotted or not. By default \code{TRUE}.}

\item{stranded}{Logical value to define whether the analyses should be performed by strand: regions in one strand will be overlapped only with regions of the same strand. The strand symbols considered are '+' and '-', any other symbol will considered in a unique separated category. Default value: \code{FALSE}.}
}
\value{
The output is the Venn Diagram in an object of class eulergram/gTree/grob/gDesc.
}
\description{
A tool to plot VennDiagrams from overlaps between bed files/regions derived from different formats. The function allows the overlap in stranded mode and can considered a specific minimal percentage of overlap between regions.
}
