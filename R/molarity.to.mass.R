#' @title Get solvent volume to make a solution with a given amount of a compound.
#'
#' @description Given a specific volume of solution wanted calculates the mass of solute necessary to obtain a certain final molarity concentration.
#'
#' @param final_concentration Numeric value for the final concentration wanted.
#' @param final_concentration_unit String to define the unit of the final concentration wanted. Available units are: "M", "mM", "uM", "nM", "pM", "fM". By default "M".
#' @param final_volume Numeric value for the final volume wanted.
#' @param final_volume_unit String to define the unit of the volume. Available units are: "L", "mL", "uL". By default "mL".
#' @param MW Numeric value for the Molecular Weigth (MW) of the compound expressed in g/mol.
#'
#' @return It returns a string with the mass of compound to use.
#'
#' @examples
#' molarity.to.mass(final_concentration = 5, final_volume = 10, MW = 215)
#'
#' @export molarity.to.mass

molarity.to.mass = function(final_concentration,
                            final_concentration_unit = "M",
                            final_volume,
                            final_volume_unit = "mL",
                            MW) { # g/mol

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


  # conversion final volume to liters
  V_unit = final_volume_unit
  V = final_volume

  if (V_unit != "L") {
    if (V_unit == "mL") {V = V * 10^-3}
    else if (V_unit == "uL") {V = V * 10^-6}
  }


  # calulate the mass
  m = MW * V * M # result in g

  # convertion of the mass
  if (m >= 1000) {mass = paste(round(m*10^-3,2), "kg", sep = "")}
  else if (0 > m & m >= 10^-3) {mass = paste(round(m*10^3,3), "mg", sep = "")}
  else if (10^-3 > m & m >= 10^-6) {mass = paste(round(m*10^6,3), "ug", sep = "")}
  else if (10^-6 > m & m >= 10^-9) {mass = paste(round(m*10^9,3), "ng", sep = "")}
  else if (10^-9 > m) {mass = paste(round(m*10^12,3), "pg", sep = "")}
  else {mass = paste(round(m,3), "g", sep = "")}

  return(mass)

} # END function
