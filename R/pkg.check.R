#' @title Check package installation.
#'
#' @description Function to check if a package is installed. It works with bioconductor or CRAN packages.
#'
#' @param package A single string indicating the name of the package to check.
#' @param archive A single string indicating the type of archive. Possible values "CRAN" and "bioconductor" (not case sensitive). Parameter without default..
#'
#' @return If the pkg is not already installed it will be installed.
#'
#' @examples
#' pkg.check("ggplot2", "cran")
#'
#' pkg.check("biomaRt", "bioconductor")
#'
#' @export pkg.check
#'
# @importFrom BiocManager install


######################
pkg.check = function(package,
                     archive) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  ### check parameters
  if (class(package) != "character" | length(package) != 1) {
    return(warning("The 'package' parameter must be a single string."))
  }


  archive = tolower(archive)
  if (class(archive) != "character" | length(archive) != 1 | !(archive %in% c("cran", "bioconductor"))) {
    return(warning("The 'archive' parameter must be a single string. Possibile values: 'cran', 'bioconductor' (not case sensitive)."))
  }



  ### install package
  if (archive == "cran") {
    # Install packages from CRAN
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      if(!require(package, character.only = TRUE)) return(warning((paste(package,"package not found."))))}

  } else {

    # Install packages from bioconductor
    if (!require(package, character.only = TRUE)) {
      BiocManager::install(package)
      if(!require(package, character.only = TRUE)) return(warning((paste(package,"package not found."))))}
  }


} # END function
