#' @title List combiner
#'
#' @description Combines two or more lists in a single one keeping the element names.
#'
#' @param list.of.lists A list of lists.
#'
#' @return It returns a list that is a combination of the lists in the input list. \cr If the list is not a nested list of list the original input is returned.
#'
#' @examples
#' combined_list = combine.lists(list.of.lists = list(list(c(1:2), c(1:3)), list("X" = c("A", "B"), "Y" = 2)))
#'
#' combined_list = combine.lists(list.of.lists = list(c(1:2), c(1:3)))
#'
#' @export combine.lists

combine.lists = function(list.of.lists) {

  #-----------------------------#
  # Check if Rseb is up-to-date #
  Rseb::actualize(update = F, verbose = F)   #
  #-----------------------------#

  # Check if the 'list.of.lists' is a list
  if (!(class(list.of.lists) == "list")) {
    return(warning("The input must be a list."))
  }


  # Check if the 'list.of.lists' is a list of list
  check = c()

  for (i in 1:length(list.of.lists)) {
    check = c(check,
              ifelse(test = class(list.of.lists[[i]]) == "list",
                     yes = T, no = F))
  }

  if (unique(check) > 1 | unique(check) == F) {
    message("The input list is not a list of list. The input list is returned as it is.")
    return(list.of.lists)
  }


  # Combine the lists
  combined.list = list()
  names = c()
  j = 0

  for (i in 1:length(list.of.lists)) {
    for (k in 1:length(list.of.lists[[i]])) {
      j = j+1
      combined.list[[j]] = list.of.lists[[i]][[k]]

      # Assign the names to the elements
      if (!is.null(names(list.of.lists[[i]])[k])) {
        names(combined.list)[j] = names(list.of.lists[[i]])[k]}
    }
  }

  return(combined.list)

} # end function
