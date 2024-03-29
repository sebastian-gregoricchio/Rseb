% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.nan_df.R
\name{is.nan_df}
\alias{is.nan_df}
\title{\code{is.nan()} applied to a data.frame}
\usage{
is.nan_df(data.frame)
}
\arguments{
\item{data.frame}{Input data.frame.}
}
\value{
It returns a matrix/array containing logic values for each element of the input data.frame. When \code{TRUE} it means that the corresponding element is a \code{NaN}.
}
\description{
Applies the function \code{is.nan()} to a full data.frame.
}
\examples{
is.nan.df(mtcars)

}
