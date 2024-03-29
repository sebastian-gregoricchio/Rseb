% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate.heterogeneity.R
\name{evaluate.heterogeneity}
\alias{evaluate.heterogeneity}
\title{Evaluate genomic heterogeneity among samples.}
\usage{
evaluate.heterogeneity(
  bigWig.list,
  peak.list,
  labels = sub(pattern = "(.*)\\\\..*$", replacement = "\\\\1", basename(bigWig.list)),
  reference.peaks = NULL,
  distribution.line.color = "#1c30a3",
  distribution.line.size = 1,
  distribution.line.type = 1,
  distribution.n.vertical.divisions = NULL,
  distribution.as.percentage = F,
  heatmap.color = "#1c30a3",
  heatmap.zMax = NA,
  heatmap.log1p.scale = TRUE,
  bar.color = distribution.line.color,
  widths.proportion = c(0.25, 1),
  heights.proportion = c(1, 1),
  min.percentage.reference = 0,
  min.percentage.test = 0,
  min.bases.overlap = 1,
  multiBigWigSummary.path = "multiBigWigSummary"
)
}
\arguments{
\item{bigWig.list}{A string vector with bigwig paths (same order than paths).}

\item{peak.list}{A list of GRanges objects (not GRglist) or data.frames or a string vector with the path to bed files (same order than bigwigs).}

\item{labels}{The labels to use for the samples (same order than bigwigs/peaks). Default: the basename of the bigWig.list.}

\item{reference.peaks}{Default: \code{NULL}, the peaks of all samples provided are merged and collapsed together.}

\item{distribution.line.color}{Color to use for the distribution line. Default: \code{"#1c30a3"} (dark blue).}

\item{distribution.line.size}{Line size of the distribution plot. Default: \code{1}.}

\item{distribution.line.type}{Line type of the distribution plot. Default: \code{1}.}

\item{distribution.n.vertical.divisions}{Number of sectors in which divide the distribution plot (vertical lines will be plotted). Default: \code{NULL} (no divisions).}

\item{distribution.as.percentage}{Logical value to define whether the distribution plot should show percentage of sample coverage rather than number of samples. Default: \code{FALSE}.}

\item{heatmap.color}{Color to use for the heatmaps; a gradient from this color to white will be used. Default: \code{"#1c30a3"} (dark blue).}

\item{heatmap.zMax}{Maximum of the heatmap scale. Default: \code{NA}.}

\item{heatmap.log1p.scale}{Logic value to define whether the heatmap scale should display log1p values. Default: \code{TRUE}.}

\item{bar.color}{Color to use for the barplot showing the fraction of reference peaks present in each sample. Default is to use the 'distribution.line.color'.}

\item{widths.proportion}{Two-elements numeric vector to be passed to \code{plot_grid} rel_width.}

\item{heights.proportion}{Two-elements numeric vector to be passed to \code{plot_grid} rel_height.}

\item{min.percentage.reference}{Numeric value within 0-100 to define which percentage of 'reference' dataset must overlap with a 'sample'. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.}

\item{min.percentage.test}{Numeric value within 0-100 to define which percentage of 'sample' datasets must overlap with a region in the 'reference' one. If the value is lower than 0 or greater than 100, will be coerced to 0 or 100 respectively. Default value: \code{0}.}

\item{min.bases.overlap}{Integer, greater than 0, value to indicate the minimal number of bases to consider as minimum overlap between two regions. Non integer values will be rounded at integer, while number lower that 1 will be coerced to 1. Default value: \code{1}.}

\item{multiBigWigSummary.path}{Path/command to run deeptools multiBigWigSummary tool. Default: \code{"multiBigWigSummary"}.}
}
\value{
The function returns a list containing:
\itemize{
  \item \code{metadata}: a list with the following elements \itemize{
                         \item \code{sample.id}: vector with the labels of the samples
                         \item \code{bigWig.file}: string vector with the path to each bigwig file
                         \item \code{bed.file}: list of peaks associated to each sample
                         \item \code{n.peaks}: vector with number of peaks in each sample
                         }
  \item \code{count.matrix}: data.frame with presence (1) or absence (0) of each peak per each sample;
  \item \code{score.matrix}: data.frame with average score at each peak per each sample;
  \item \code{plot.list}: list of single separated plots: counts.distribution, fraction.peaks.per.sample, scores.heatmap;
  \item \code{multiplot}: the multiplot generated from the plot.list.
 }
}
\description{
This tools evaluates what is the fraction of peaks covered by each sample provided in a peaks dataset obtained by merging all the peaks together or provided by the user. The peaks in the reference dataset are ranked by number of samples in which are present and average score all over the samples. This function uses the deeptools function multiBigwigSummary.
}
\details{
To know more about the deepTools's function \code{multiBigwigSummary} see the package manual at the following link: \cr \url{https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html}.
}
