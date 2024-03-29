% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic.tracks.R
\name{genomic.tracks}
\alias{genomic.tracks}
\title{Genomic tracks plotter}
\usage{
genomic.tracks(
  tracks,
  genomic.region,
  genome,
  track.labels = NULL,
  track.labels.fontzise = 5,
  track.labels.position = c(-0.1, 0),
  track.colors = "#000000",
  grouping = NULL,
  gene.annotation.color = "darkblue",
  expand.bed = TRUE,
  arcs.direction = "down",
  fraction.arc.base = 0.025,
  highlight.bed = NULL,
  highlight.color = "yellow",
  highlight.transparency = 0.15,
  missing.data.as.zero.bw = FALSE,
  smooth.bigWig.signal = TRUE,
  smooth.bigWig.loess.span = 0.05,
  plot.bigWig.area = TRUE,
  bigWig.range.label.size = 2.5,
  score.bed.shadow = FALSE,
  height.ratios = NULL,
  width.ratios = c(1, 5)
)
}
\arguments{
\item{tracks}{A vector indicating the list of full paths of the files/tracks/signals to plot. Supported formats: bed/bd/narrowPeak/broadPeak, bw/bigWig/bigwig, bedpe.}

\item{genomic.region}{An atomic string indicating the genomic region into which restrict the final plot in the format 'chr1:1234-5678'.}

\item{genome}{An atomic string indicating the genome to use for the annotations. Allowed values are:
\itemize{
  \item \code{hg19}: loads an 'EnsDb' object from the library \code{EnsDb.Hsapiens.v75};
  \item \code{hg38}: loads an 'EnsDb' object from the library \code{EnsDb.Hsapiens.v86};
  \item \code{mm10}: loads an 'EnsDb' object from the library \code{EnsDb.Mmusculus.v79};
  \item \emph{custom 'EnsDb' object}: provide an 'EnsDb' object manually generated; visit the page \url{https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#102_building_annotation_packages} for more information.
 }}

\item{track.labels}{A vector indicating the labels to use for each track (genome annotation track excluded). By default \code{NULL}: the file base-name will be used.}

\item{track.labels.fontzise}{A numerical value to indicate the font size of the track labels. Default value \code{5}.}

\item{track.labels.position}{A two-element numeric vector passed to \code{xlim} function for the the definition of the frame size of the track labels. Default value \code{c(-0.1, 0)}.}

\item{track.colors}{A string vector indicating the color to use for each track (genome annotation track excluded). If only one value is provided it will be used for all the tracks. Default value \code{"#000000"} ("black").}

\item{grouping}{A single numerical vector or a list of numeric vectors. Each list's element indicates the indexes corresponding to the tracks (1 = first track, 2 = second track, etc) for which the y-axes should be normalized. Each element will be taken into account in the order. Default value \code{NULL}.}

\item{gene.annotation.color}{A string indicating the color to use for the genome annotation track.}

\item{expand.bed}{A logical value to define whether overlapping regions in a bed should be plotted on different levels. Default \code{TRUE}.}

\item{arcs.direction}{A string indicating the direction on which arcs should be plotted for bedpe files. Available options \code{"up"} or \code{"down"}. Default value \code{"down"}.}

\item{fraction.arc.base}{A numerical value indicating the fraction of total plot height to be used as arc base thickness. By default \code{0.025} (2.5\% of the track height).}

\item{highlight.bed}{Either a string indicating the full path to a bed file or a data.frame in BED3 format (chr, start, end) containing regions that should be highlighted in the plot. Regions included in the genomic range will be automatically selected. By default \code{NULL}.}

\item{highlight.color}{A string indicating the color to use for the regions to highlight in the plot. By default \code{'yellow'}.}

\item{highlight.transparency}{A numerical value indicating the transparency (alpha) to use for the highlighted regions. Default value \code{0.15}.}

\item{missing.data.as.zero.bw}{A logical value to define wthere missing data in the bigWigs should be converted to zeros. Default \code{FALSE}.}

\item{smooth.bigWig.signal}{Logical value to indicate whether the bigWig signals should be smoothed (by loess x ~ y function). By default \code{TRUE}.}

\item{smooth.bigWig.loess.span}{Numerical value to indicate the span value for the loess function used to smooth bigWig signals. By default \code{0.05}.}

\item{plot.bigWig.area}{Logical value to indicate whether the bigWig profile should be filled or not. If \code{FALSE} only the signal outline will be plotted. By default \code{TRUE}.}

\item{bigWig.range.label.size}{A numerical value to indicate the font size of the bigWig signal range. Default value \code{2.5}.}

\item{score.bed.shadow}{Logical value to define whether the filling intensity of the bed segments should reflect the score of each signal. By default \code{FALSE}.}

\item{height.ratios}{Numerical vector of relative track heights, passed to 'rel_heights' parameter of \code{cowplot::plot_grid()}. For example, in a two-row grid, rel_heights = c(2, 1) would make the first column twice as wide as the second column. Value \code{1} indicates that all the tracks should have the same size. By default \code{NULL}, automatic ratios will be computed by this function.}

\item{width.ratios}{Numerical vector of relative labels vs tracks widths, passed to 'rel_widths' parameter of \code{cowplot::plot_grid()}. For example, in a two-column grid, rel_widths = c(2, 1) would make the first column twice as wide as the second column. Value \code{1} indicates that all the tracks should have the same size. By default \code{c(1,5)} (1 label : 5 tracks).}
}
\value{
The function returns a named list containing:
\itemize{
  \item \code{configuration}: data.frame with the parameters used to build the plot(s);
  \item \code{highlighted.region}: data.frame with the regions used for the highlighting;
  \item \code{single.track.list}: a named list containing each single track plot used for the creation of the multi.track.plot;
  \item \code{single.label.plot.list}: a named list containing each single track label plot used for the creation of the multi.track.plot;
  \item \code{multi.track.plot}: the assembled multi.track labelled plot.
 }
}
\description{
The functions allows to plot different types of genomic data (bigWig, bed, bedpe) at a specific genomic region. It is possible to highlight specific regions and the gene annotations are plotted automatically at the bottom of all the tracks.
}
