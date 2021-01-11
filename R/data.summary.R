#' @title Statistical data summary generator
#'
#' @description Produces a table with a summary of the statistics for a specific column of an input data.frame by a group of values defined by a group defined by another column.
#'
#' @param data Input data.frame to be analyzed.
#' @param variable A string with the name of the column to be analyzed.
#' @param group.names A string with the name of the column indicating the groups.
#'
#' @return It returns a list that is a combination of the lists in the input list. \cr If the list is not a nested list of list the original input is returned.
#'
#' @examples
#' data.summary(data = mtcars, variable = "mpg", group.names = "disp")
#'
#' @export data.summary
#'
# @importFrom plyr ddply


data.summary = function(data,
                        variable,
                        group.names){

  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x)),
      median = median(x[[col]], na.rm=TRUE))
  }

  data_sum = plyr::ddply(data,
                         group.names,
                         .fun=summary_func,
                         variable)
  #data_sum = rename(data_sum, c("mean" = variable))
  return(data_sum)
}
