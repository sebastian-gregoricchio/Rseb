#' @title function to automatically update the R packages.
#'
#' @description Automatically updates the R packages from CRAN and BioConductor repositories.
#'
#' @param ask Logical indicating whether to ask the user to select packages before they are downloaded and installed, or the character string \code{"graphics"}, which brings up a widget to allow the user to (de-)select from the list of packages which could be updated. (The latter value only works on systems with a GUI version of \code{select.list}, and is otherwise equivalent to \code{ask = TRUE}). By default \code{FALSE}.
#'
#' @return Nothing. The packages will be updated.
#'
#' @examples
#' update_pkgs()
#'
#' @export update_pkgs

update_pkgs = function(ask = FALSE) {

  # Update CRAN libraries
  update.packages(ask = ask)


  # Update BIOCONDUCTOR libraries
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(ask = ask)

} # end function
