#' @title Get session info and package versions.
#'
#' @description Retrieves the information of the current session and the version of the packages loaded.
#'
#' @param return.session Logic value to define if to save the session info. By default \code{FALSE}.
#' @param print.versions Logic value to define if to print the session and version info. By default \code{TRUE}.
#' @param return.versions Logic value to define if to save package versions info. By default \code{FALSE}.
#'
#' @return If \code{return.session} and/or \code{return.versions} \code{TRUE} a list with these informations is returned. Otherwise nothing is returned.
#'
#' @export pkg.version

pkg.version = function(return.session = F,
                       print.versions = T,
                       return.versions = F) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F,   #
                  verbose = F)  #
  #-----------------------------#

  # Check if Rseb is attached, if not, attach it. then it will be eventually removed again
  was.attached = "Rseb" %in% (.packages())
  require("Rseb", quietly = T)

  session = sessionInfo()

  list.v = list(`R version:` = session$R.version$version.string,
                `Base packages:` = sort(session$basePkgs),
                `Other attached packages:` = sort(paste(names(session$otherPkgs),
                                                        installed.packages()[names(session$otherPkgs), "Version"],
                                                        sep = "_")),
                `Loaded via a namespace (and not attached):` = sort(paste(names(session$loadedOnly),
                                                                          installed.packages()[names(session$loadedOnly), "Version"],
                                                                          sep = "_")))
  # Indicating when a section is empty
  for (i in 1:length(list.v)) {
    if (length(list.v[[i]]) == 0) {
      (list.v[[i]] = "No packages loaded in this section")}
  }

  versions = list(session = session,
                  versions = list.v)

  if (print.versions == T) {print(versions$versions, quote = F)}

  if (return.session == T | return.versions == T) {
    if (return.session == T & return.versions == F) {
      return(versions$session)} else if (return.session == F & return.versions == T) {
        return(versions$versions)} else {return(versions)}
  }

  # Detach Rseb if it was not already attached previously
  if (was.attached == F) {detach("package:Rseb", unload = TRUE)}

} # end function
