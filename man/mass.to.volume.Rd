% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mass.to.volume.R
\name{mass.to.volume}
\alias{mass.to.volume}
\title{Get solvent volume to make a solution with a given amount of a compound.}
\usage{
mass.to.volume(
  final_concentration,
  final_concentration_unit = "M",
  mass,
  mass_unit = "g",
  MW
)
}
\arguments{
\item{final_concentration}{Numeric value for the final concentration wanted.}

\item{final_concentration_unit}{String to define the unit of the final concentration wanted. Available units are: "M", "mM", "uM", "nM", "pM", "fM". By default "M".}

\item{mass}{Numeric value for the solute mass ammount.}

\item{mass_unit}{String to define the unit of the mass. Available units are: "kg", "g", "mg", "ug", "ng". By default "g".}

\item{MW}{Numeric value for the Molecular Weigth (MW) of the compound expressed in g/mol.}
}
\value{
It returns a string with the volume of solvent to use.
}
\description{
Given a specific ammount of solute calculates the volume of solvent necessary to obtain a certain final molarity concentration.
}
\examples{
mass.to.volume(final_concentration = 5, mass = 10, MW = 215)

}
