% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restore_packages.R
\name{restore_packages}
\alias{restore_packages}
\title{Restores packages installed from a .rda file.}
\usage{
restore_packages(rda_file)
}
\arguments{
\item{rda_file}{Path to the .rda from which get the information for the packages to re-install.}
}
\value{
If it was not possible to re-install al packages, the list of not restored packages will be returned.
}
\description{
Installs the packages contained in a .rda file. This file can be generated by the \code{\link{store_packages}} function of this package.
}
