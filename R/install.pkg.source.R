#' @title Package installer from source archive.
#'
#' @description Allows the installation of R packages using the source archive file.
#'
#' @param pkg.path String to define the path for the archive file to be installed.
#'
#' @return No returned value. The package required will be installed.
#'
#' @export install.pkg.source

install.pkg.source = function(pkg.path) {
  install.packages(pkg.path, repos = NULL, type="source")
}
