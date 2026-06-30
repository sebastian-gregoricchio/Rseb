.get_option <- function(name, default) {
  getOption(paste0("Rseb.", name), default)
}


#' @title Package options
#'
#' @description The \code{Rseb} package uses the following global options
#' (set via \code{\link[base:options]{options()}}):
#'
#' \describe{
#'   \item{\code{Rseb.update_check}}{
#'     Logical. Whether Rseb checks for package updates.
#'     Default is \code{TRUE}.
#'   }
#' }
#'
#' @name Rseb-options
#' @docType data
#' @keywords internal
NULL