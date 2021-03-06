% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkg.version.R
\name{pkg.version}
\alias{pkg.version}
\title{Get session info and package versions.}
\usage{
pkg.version(return.session = F, print.versions = T, return.versions = F)
}
\arguments{
\item{return.session}{Logic value to define if to save the session info. By default \code{FALSE}.}

\item{print.versions}{Logic value to define if to print the session and version info. By default \code{TRUE}.}

\item{return.versions}{Logic value to define if to save package versions info. By default \code{FALSE}.}
}
\value{
If \code{return.session} and/or \code{return.versions} \code{TRUE} a list with these informations is returned. Otherwise nothing is returned.
}
\description{
Retrieves the information of the current session and the version of the packages loaded.
}
