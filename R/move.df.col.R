#' @title Function to change easily the order of specific columns in a data.frame.
#'
#' @description Allows to change the position of a column in a data.frame using other columns as reference.
#'
#' @param data.frame An input data.frame.
#' @param move.command A string containing the moving command. The command is formed as follows: "columnA movingCommand columnB". The basic options are: "first", "last", "before", "after". Compounded moves must be separated by a semicolon. Example: \code{"g first; a last; e before c"}.
#'
#' @return It returns the original data.frame but with the columns moved as demanded.
#'
#' @examples
#' new.mtcars = move.df.col(mtcars, "mpg last")
#'
#' new.mtcars = move.df.col(mtcars, "wt before carb")
#'
#' new.mtcars = move.df.col(mtcars, "am before carb; cyl first")
#'
#' @references \url{https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping}
#'
#' @export move.df.col

move.df.col = function (data.frame,
                         move.command) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  move.command = lapply(strsplit(strsplit(move.command, ";")[[1]],
                                 ",|\\s+"), function(x) x[x != ""])

  movelist = lapply(move.command, function(x) {
    Where = x[which(x %in% c("before", "after", "first",
                              "last")):length(x)]
    ToMove = setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec = names(data.frame)
  for (i in seq_along(movelist)) {
    temp = setdiff(myVec, movelist[[i]][[1]])
    A = movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba = movelist[[i]][[2]][2]
      if (A == "before") {
        after = match(ba, temp) - 1
      }
      else if (A == "after") {
        after = match(ba, temp)
      }
    }
    else if (A == "first") {
      after = 0
    }
    else if (A == "last") {
      after = length(myVec)
    }
    myVec = append(temp, values = movelist[[i]][[1]], after = after)
  }

  return(data.frame[myVec])
}
