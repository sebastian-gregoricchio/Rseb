#' @title \code{Rseb} updates verification
#'
#' @description It verifies if \code{Rseb} is up-to-date and installs it when required.
#'
#' @param update Logical value to define whether update the \code{Rseb} package. By default \code{TRUE}.
#' @param verbose Logical value to define whether print messages. By default \code{TRUE}.
#' @param force Logical value to define whether to force the installation of \code{Rseb} even though already up-to-date. Parameter passed to \code{devtools::install_github()}. By default \code{FALSE}.
#'
#' @return Warnings and/or messages. Installation of the latest version of \code{Rseb} if required.
#'
#' @details This function will check for internet availability.
#'
#' @export actualize


actualize = function(update = TRUE,
                     verbose = TRUE,
                     force = FALSE) {

  # Check internet connection
  if (isFALSE(curl::has_internet())) {return(invisible(NULL))}

  # Store Rseb package versions
  check = suppressMessages(rvcheck::check_github("sebastian-gregoricchio/Rseb"))
        # When up-to-date == NA it means that it is a developmental version higher then the last release
        if (is.na(check$up_to_date)) {return(invisible(NULL))}

  # Check if Rseb is up-to-date
  if (check$up_to_date == FALSE) {
      return(warning(paste("The Rseb package is not up-to-date. Installed version",
                           check$installed_version, "-->", check$latest_version,
                           "available. \nTo update Rseb type: Rseb::actualize()")))
  }

  # Update Rseb if required
  if (check$up_to_date == FALSE & update == TRUE) {
    devtools::install_github("sebastian-gregoricchio/Rseb", force = force)
  } else if (verbose == T) {
    return(message(paste("Rseb's latest version (v", check$installed_version, ") is already installed.", sep = "")))}

} # END function
