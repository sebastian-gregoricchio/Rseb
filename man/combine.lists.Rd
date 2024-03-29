% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.lists.R
\name{combine.lists}
\alias{combine.lists}
\title{List combiner}
\usage{
combine.lists(list.of.lists)
}
\arguments{
\item{list.of.lists}{A list of lists.}
}
\value{
It returns a list that is a combination of the lists in the input list. \cr If the list is not a nested list of list the original input is returned.
}
\description{
Combines two or more lists in a single one keeping the element names.
}
\examples{
combined_list = combine.lists(list.of.lists = list(list(c(1:2), c(1:3)), list("X" = c("A", "B"), "Y" = 2)))

combined_list = combine.lists(list.of.lists = list(c(1:2), c(1:3)))

}
