#' @title Get solvent volume to make a solution with a given amount of a compound.
#'
#' @description Given a specific ammount of solute calculates the volume of solvent necessary to obtain a certain final molarity concentration.
#'
#' @param final_concentration Numeric value for the final concentration wanted.
#' @param final_concentration_unit String to define the unit of the final concentration wanted. Available units are: "M", "mM", "uM", "nM", "pM", "fM". By default "M".
#' @param mass Numeric value for the solute mass ammount.
#' @param mass_unit String to define the unit of the mass. Available units are: "kg", "g", "mg", "ug", "ng". By default "g".
#' @param MW Numeric value for the Molecular Weigth (MW) of the compound expressed in g/mol.
#'
#' @return It returns a string with the volume of solvent to use.
#'
#' @examples
#' mass.to.volume(final_concentration = 5, mass = 10, MW = 215)
#'
#' @export mass.to.volume

mass.to.volume = function(final_concentration,
                          final_concentration_unit = "M",
                          mass,
                          mass_unit = "g",
                          MW) { # g/mol

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # convertion final concentration to molar
  M_unit = final_concentration_unit
  M = final_concentration

  if (M_unit != "M") {
    if (M_unit == "mM") {M = M * 10^-3}
    else if (M_unit == "uM") {M = M * 10^-6}
    else if (M_unit == "nM") {M = M * 10^-9}
    else if (M_unit == "pM") {M = M * 10^-12}
    else if (M_unit == "fM") {M = M * 10^-15}
  }

  # convertion of mass to g
  m = mass

  if (mass_unit != "g") {
    if (mass_unit == "kg") {m = m * 10^3}
    else if (mass_unit == "mg") {m = m * 10^-3}
    else if (mass_unit == "ug") {m = m * 10^-6}
    else if (mass_unit == "ng") {m = m * 10^-9}
  }

  # calulate the mass
  V = m / (MW * M) # result in L


  # calculate the final volume
  if (0.5 > V & V >= 10^-3) {volume = paste(round(V*10^3,2), "mL", sep = "")}
  else if (10^-3 > V) {volume = paste(round(V*10^6,1), "uL", sep = "")}
  else {volume = paste(round(V,2), "L", sep = "")}


  return(volume)

} # END function
