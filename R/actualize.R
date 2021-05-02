#' @title \code{Rseb} updates verification
#'
#' @description It verifies if \code{Rseb} is up-to-date and installs it when required.
#'
#' @param update Logical value to define whether update the \code{Rseb} package. By default \code{TRUE}.
#' @param verbose Logical value to define whether print messages. By default \code{TRUE}.
#' @param force Logical value to define whether to force the installation of \code{Rseb} even though already up-to-date. Parameter passed to \code{devtools::install_github()}. By default \code{FALSE}.
#' @param build.manual Logic value to define whether to build the manual. By default \code{TRUE}.
#' @param build.vignettes Logic value to define whether to build the vignettes. By default \code{TRUE}.
#'
#' @return Warnings and/or messages. Installation of the latest version of \code{Rseb} if required.
#'
#' @details This function will check for internet availability.
#'
#' @export actualize


actualize = function(update = TRUE,
                     verbose = TRUE,
                     force = FALSE,
                     build.manual = TRUE,
                     build.vignettes = TRUE) {

  # Check internet connection
  if (isFALSE(curl::has_internet())) {return(invisible(NULL))}

  # Store Rseb package versions
  check = suppressMessages(rvcheck::check_github("sebastian-gregoricchio/Rseb"))
  # When up-to-date == NA it means that it is a developmental version higher then the last release
  if (is.na(check$up_to_date)) {return(invisible(NULL))}

  # Check if Rseb is up-to-date and print a message if not up-to-date
  if (check$up_to_date == F & update == F) {
    message = paste("| The 'Rseb' package is not up-to-date. Installed version v",
                    check$installed_version, " --> v", check$latest_version,
                    " available. |", sep="")
    command = "| To update Rseb type: Rseb::actualize()"

    message(paste("\033[0;37;44m",
                  paste(rep("-", nchar(message)), collapse = ""),
                  "\n", message, "\n",
                  command, paste(rep(" ", (nchar(message) - nchar(command)) - 1),  collapse = ""), "|",
                  "\n", paste(rep("-", nchar(message)), collapse = ""),
                  "\033[0m",
                  sep = ""))
    return(invisible(NULL))
  }


  # Update Rseb if required
  if (check$up_to_date == F & update == T) {
    message(paste("\033[0;37;44m The 'Rseb' package is not up-to-date.\n Rseb will be updated to the last release (v",
                  check$installed_version, " --> v", check$latest_version, "). \033[0m", sep=""))
    devtools::install_github("sebastian-gregoricchio/Rseb", build_manual = build.manual, build_vignettes = build.vignettes)}
  else if (check$up_to_date == T) {
    if (force == T) {devtools::install_github("sebastian-gregoricchio/Rseb", force = force, build_manual = build.manual, build_vignettes = build.vignettes)}
    else if (verbose == T) {return(message(paste("Rseb's latest version (v", check$installed_version, ") is already installed.", sep = "")))}
  }


} # END function
