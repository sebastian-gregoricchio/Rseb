#' @title Grep a pattern in a full data.frame.
#'
#' @description The function helps to define which rows of an input data.frame contain a specific patter.
#'
#' @param data.frame Input data.frame.
#' @param pattern Character string containing a regular expression (or character string for \code{fixed = TRUE}) to be matched in the given character vector. Coerced by \code{as.character} to a character string if possible. If a character vector of length 2 or more is supplied, the first element is used with a warning. Missing values are allowed except for \code{regexpr} and \code{gregexpr}.
#' @param ignore.case If \code{FALSE}, the pattern matching is case sensitive and if \code{TRUE}, case is ignored during matching. By default \code{FALSE}.
#' @param perl Logical value to define if Perl-compatible regexps should be used. By default \code{FALSE}.
#' @param fixed Logical value to define if the pattern is a string to be matched as is. Overrides all conflicting arguments. By default \code{FALSE}.
#' @param useBytes Logical value to define if the matching is done byte-by-byte rather than character-by-character. By default \code{FALSE}.
#'
#' @return It will be return a logic vector with an element per each row of the data.frame. The value is \code{TRUE} when the patter is found at least once in the corresponding data.frame row.
#'
#' @examples
#' iris = iris %>% filter(grepl.data.frame(iris, pattern = "setosa"))
#'
#' @export grepl.data.frame


grepl.data.frame =
  function(data.frame,
           pattern,
           ignore.case = FALSE,
           perl = FALSE,
           fixed = FALSE,
           useBytes = FALSE) {

  df = data.frame(data.frame, stringsAsFactors = F)

  for (i in 1:length(colnames(df))) {
    df[,i] = grepl(x = df[,i],
                   pattern = pattern,
                   ignore.case = ignore.case,
                   perl = perl,
                   fixed = fixed,
                   useBytes = useBytes)
  }

  check = c()
  names(df) = NULL
  for (i in 1:nrow(df)) {
    if (length(unique(unlist(df[i,]))) == 2) {
      values = unique(unlist(df[i,]))
      check[i] = values[1] | values[2]
    } else (check[i] = unique(unlist(df[i,])))
  }

  return(check)
}
