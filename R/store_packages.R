#' @title Stores the information of installed packages in a .rda file.
#'
#' @description Saves the list of all the installed packages in a .rda file. This file can be used to restore the packages from a computer to another or after installation of a new R version by the function \code{\link{restore_packages}} of this package.
#'
#' @param output_directory Path to the directory in which export the .rda file. By default \code{<working.directory>}.
#'
#' @return Nothing is returned. An .rda file will be exported at the \code{output_directory} indicated.
#'
#' @export store_packages

store_packages = function(output_directory = getwd()) {
  # stores a list of your currently installed packages
  tmp = installed.packages()
  installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])

  date = Sys.Date()
  file_name = paste(c(output_directory,"/installed_packages_",as.character(as.data.frame(Sys.info())["user",]),"_",as.character(Sys.Date()),".rda"),
                    sep = "", collapse = "")

  save(installedpackages, file=file_name)
}
