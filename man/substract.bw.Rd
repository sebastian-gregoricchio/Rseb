% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/substract.bw.R
\name{substract.bw}
\alias{substract.bw}
\title{Combination of two or more list in a unique one.}
\usage{
substract.bw(bw1, bw2, return.substracted.bw = T, substracted.bw.file = NULL)
}
\arguments{
\item{bw1}{Full path to the first bigWig (the second one will be substracted to this one).}

\item{bw2}{Full path to the second bigWig (it will be substracted to the first one).}

\item{return.substracted.bw}{Logic value to define whether return the resulting bigWig as GRanges object. By default \code{TRUE}.}

\item{substracted.bw.file}{String for the path of the resulting bigwig file to be exported. \cr By default \code{NULL}, any file will be exported.}
}
\value{
If required a subtraction bigWig is returned as GRanges object. The resulting bigWig can be also directly exported.
}
\description{
Combines two or more lists in a single one keeping the elements names
}
