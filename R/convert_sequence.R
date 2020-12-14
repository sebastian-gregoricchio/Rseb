#' @title Nucleic acid sequences converter.
#'
#' @description Obtains de complementary, reverse complementary or the reverse of a DNA/RNA sequence.
#'
#' @param sequence A string containing the sequence to be converted. By default \code{NULL}, it returns an help for the mode.
#' @param mode A string value to define the modality of convertion. Possible options: \cr - Reverse complement = revComp | RC | rc | reverseComplement \cr - Reverse            = rev     | R  | r  | reverse \cr - Complement         = comp    | C  | c  | complement. \cr By default \code{"not specified"}, it returns an help for the mode.
#' @param nucleic.acid A string to define the type of nucleic acid to which the input sequence belongs. Available options "DNA", default value, or "RNA".
#'
#' @return It returns a string with the converted sequence.
#'
#' @examples
#' convert_sequence(sequence = "AATTTCCCGTCGAT",
#'                  mode = "reverse",
#'                  nucleic.acid = "DNA")
#'
#' @export convert_sequence
#'
#' @import Biostrings

convert_sequence = function(sequence = NULL,
                            mode = "not specified",
                            nucleic.acid = "DNA") {

  ###### Install required packages  ######
  pkg = "Biostrings"
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    if(!require(pkg, character.only = TRUE)) stop(paste(pkg,"package not found."))}


  ###### Parameters check ######
  if (class(sequence) != "character" | is.null(sequence)) {
    return(warning("The sequence must be a string of class <character>"))
  }

  if (mode == "not specified" | !(mode %in% c("revComp", "RC", "rc", "reverseComplement",
                                              "rev", "R", "r", "reverse",
                                              "comp", "C", "c", "complement")) |
      class(mode) != "character") {
    return(warning(c("\n Transformation mode is not available or not specified. \n \n",
                     "Possible choices: \n",
                     "- Reverse complement = revComp | RC | rc | reverseComplement \n",
                     "- Reverse            = rev     | R  | r  | reverse \n",
                     "- Complement         = comp    | C  | c  | complement")))
    }

  if (!(nucleic.acid %in% c("DNA", "dna", "RNA", "rna")) | class(nucleic.acid) != "character") {
    return(warning(c("Nucleic acid type not recognized. \n \n",
                     "Possible choices: DNA, dna | RNA, rna")))
  }

  ######################################################

  ###### Check library ######

  if (!("Biostrings" %in% installed.packages())) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings")
  }

  require("Biostrings")

  ######################################################


  ######## Sequence transformation
  if (nucleic.acid == "DNA") {
    acid.sequence = DNAString(toupper(sequence))} else {
      acid.sequence = RNAString(toupper(sequence))}

  if (mode %in% c("revComp", "RC", "rc", "reverseComplement")) {
    return(toString(reverseComplement(acid.sequence)))}

  if (mode %in% c("rev", "R", "r", "reverse")) {
    return(toString(reverse(acid.sequence)))}

  if (mode %in% c("comp", "C", "c", "complement")) {
    return(toString(complement(acid.sequence)))}

} #end
