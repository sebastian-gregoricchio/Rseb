#' @title Intersect two or more bed files (by \code{bedtools intersect} function).
#'
#'
#' @description This function runs a command line that uses \code{bedtools intersect} to intersect one or more .bed files.
#'
#'
#' @param a A single string defining the BAM/BED/GFF/VCF file “A”. Each feature in A is compared to B in search of overlaps. Use “stdin” if passing A with a UNIX pipe.
#' @param b A character vector with one or more BAM/BED/GFF/VCF file(s) “B”. It could be also a single string containing wildcard (*) character(s).
#' @param outputFileName Full path to output file name. By default \code{<working.directory>/intersected.bed}.
#'
#' @param abam Logic value to define if file A is a BAM. Each BAM alignment in A is compared to B in search of overlaps. By default \code{FALSE}.
#' @param ubam Logic value to define if to write the output as uncompressed BAM. The default is to write compressed BAM output (\code{ubam = FALSE}).
#' @param bed Logic value to define whether to write output as BED when using a BAM input \code{abam = TRUE}. The default is to write output in BAM (\code{bed = FALSE}).
#'
#' @param wa Logic value to define if to write the original entry in A for each overlap. By default \code{FALSE}.
#' @param wb Logic value to define if to write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r. By default \code{FALSE}.
#' @param loj Logic value to define if to perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B. By default \code{FALSE}.
#' @param wo Logic value to define if to write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r. By default \code{FALSE}.
#' @param wao Logic value to define if to write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r. By default \code{FALSE}.
#' @param u Logic value to define if to write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B. Restricted by -f and -r. By default \code{FALSE}.
#' @param c Logic value to define if to for each entry in A, report the number of hits in B while restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted -f, -F, -r, and -s. By default \code{FALSE}.
#' @param C Logic value to define if to for each entry in A, separately report the number of overlaps with each B file on a distinct line. Reports 0 for A entries that have no overlap with B. Overlaps restricted by -f, -F, -r, and -s. By default \code{FALSE}.
#' @param v Logic value to define if to only report those entries in A that have no overlap in B. Restricted by -f and -r.
#'
#' @param f Numeric value defining the minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp). By default \code{NULL}.
#' @param F. Numeric value defining the minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp). By default \code{NULL}.
#' @param r Numeric value defining the required reciprocal fraction of overlap for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90\% of A and that A also overlaps at least 90\% of B. By default \code{NULL}.
#' @param e Numeric value defining the minimum fraction to be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90\% of A is covered OR 10\% of B is covered. Without -e, both fractions would have to be satisfied. By default \code{NULL}.
#' @param s Logic value to define if to force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand. By default \code{FALSE}.
#' @param S Logic value to define if to require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand. By default \code{FALSE}.
#'
#' @param split Logic value to define if to treat “split” BAM (i.e., having an “N” CIGAR operation) or BED12 entries as distinct BED intervals. By default \code{FALSE}.
#' @param sorted Logic value to define, for very large B files, if to invoke a “sweeping” algorithm that requires position-sorted input. When using -sorted, memory usage remains low even for very large files. By default \code{FALSE}. It is possible to sort a bed file on terminal by \code{(sort -k1,1 -k2,2n unsorted.bed > sorted.bed)} or by the function \link{sort.bed}.
#' @param g Specify a genome file the defines the expected chromosome order in the input files for use with the -sorted option. By default \code{NULL}.
#'
#' @param srun Logic value to define whether the command should be run in \code{srun} mode. By default \code{FALSE}.
#' @param intersect.bedtools.command String to define the command to use to recall the \code{bedtools intersect} function. An example: "/home/user/anaconda3/bin/intersectBed". By default \code{"intersectBed"}.
#'
#' @param return.command Logic value to define whether to return the string corresponding to the command for bedtools. By default \code{FALSE}.
#' @param return.bed Logic value to define whether to return the resulting bed as data.frame. By default \code{FALSE}. Parameter not active when inputs are bam files.
#' @param delete.output Logic value to define whether to delete the exported intersected bed file. By default \code{FALSE}. Parameter active only when \code{return.bed = TRUE}. Useful when is sufficient to get the result as a data.frame without saving it.
#' @param run.command Logic value to define whether to run the the command line on system terminal and generate the bed resulting from the intersection. By default \code{TRUE}.
#'
#'
#' @return The function generates the files indicated by the output parameters. If required the command line used and/or the resulting intersected bed file. If both outputs are required, the output will be a named list with two values: "command" and "intersected.bed".
#'
#' @examples
#' intersect.bedtools(a = bed_file1.bed,
#'                    b = c("bed_file2.bed", "bed_file3.bed"),
#'                    wb = TRUE,
#'                    intersect.bedtools.command = "/home/user/anaconda3/bin/intersectBed")
#'
#' intersect.bedtools(a = bed_file1.bed,
#'                    b = c("bed_file2.bed", "bed_file3.bed"),
#'                    wa = TRUE,
#'                    return.bed = TRUE,
#'                    delete.output = T,
#'                    intersect.bedtools.command = "/home/user/anaconda3/bin/intersectBed")
#'
#' @details To know more about the \code{bedtools intersect} function see the package manual at the following link: \cr \url{https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}.
#'
#' @export intersect.bedtools
#'
#' @import data.table fread



intersect.bedtools =
  function(
    # basic input/output parameters
    a,
    b,
    outputFileName = paste(getwd(), "intersected.bed", sep = "/"),

    # options for BAM files
    abam = FALSE,
    ubam = FALSE,
    bed = FALSE,

    # resulting intersection options
    wa = FALSE,
    wb = FALSE,
    loj = FALSE,
    wo = FALSE,
    wao = FALSE,
    u = FALSE,
    c = FALSE,
    C = FALSE,
    v = FALSE,

    # options for overlaps restriction
    f = NULL,
    F. = NULL,
    r = NULL,
    e = NULL,
    s = FALSE,
    S = FALSE,

    # files handling options
    split = FALSE,
    sorted = FALSE,
    g = NULL,

    # path to bedtools
    srun = FALSE,
    intersect.bedtools.command = "intersectBed",

    # others parameters for this function in R
    return.command = FALSE,
    return.bed = FALSE,
    delete.output = FALSE,
    run.command = TRUE
  ) {


    ######################################################################################
    # Create function to add single quote to string variables
    add.quotes = function(x) {return(sapply(x, function(x){paste("'", x, "'", sep = "")}, USE.NAMES = F))}


    ###### Generation of the string to run the command on bash
    if (srun == TRUE) {
      command = paste("srun", intersect.bedtools.command)
    } else {command = intersect.bedtools.command}


    ###### Check 'a' and 'b' input files
    if (length(a) != 1) {return(warning("Only one A input bed is allowed."))}
    if (class(a) != "character" | class(b) != "character") {return(warning("The A and B bed inputs must be strings, or a string vector for B input, with the full path to the input files."))}


    ###### BAM parameters
    if (abam == T) {
      command = paste(command,
                      "-abam", add.quotes(a),
                      "-b", paste(add.quotes(b), collapse = " "))
      if (ubam == T & bed == F) {
        command = paste(command, "-ubam")
      } else if (ubam == F & bed == T) {
        command = paste(command, "-bed")
        } else if (ubam == T & bed == T) {
          return(warning("The parameters 'ubam' and 'bed' cannot be TRUE simultaneously."))
        }
    } else { # when abam == FALSE
      command = paste(command,
                      "-a", add.quotes(a),
                      "-b", paste(add.quotes(b), collapse = " "))
    }


    ###### Resulting intersection options
    if (wa == T) {command = paste(command, "-wa")}
    if (wb == T) {command = paste(command, "-wb")}
    if (loj == T) {command = paste(command, "-loj")}
    if (wo == T) {command = paste(command, "-wo")}
    if (wao == T) {command = paste(command, "-wao")}
    if (u == T) {command = paste(command, "-u")}
    if (c == T) {command = paste(command, "-c")}
    if (C == T) {command = paste(command, "-C")}
    if (v == T) {command = paste(command, "-v")}


    ###### options for overlaps restriction
    if (!is.null(f)) {command = paste(command, f)}
    if (!is.null(F.)) {command = paste(command, F.)}
    if (!is.null(r)) {command = paste(command, r)}
    if (!is.null(e)) {command = paste(command, e)}
    if (s == T) {command = paste(command, "-s")}
    if (S == T) {command = paste(command, "-S")}


    ###### files handling options
    if (split == T) {command = paste(command, "-split")}
    if (sorted == T) {command = paste(command, "-sorted")}
    if (!is.null(g) & length(g) == 1) {command = paste(command, add.quotes(g))}


    ###### Add output command
    command = paste(command, ">", add.quotes(outputFileName))


    ###### Run and/or export the command
    if (run.command == T) {
      system(command)
      message("The intersection has been performed.")}

    if (return.command == T & return.bed == F) {
      return(command)
    } else if (return.command == T & return.bed == T) {
        bed = as.data.frame(data.table::fread(outputFileName))
        if (delete.output == T) {file.remove(outputFileName)}
        return(list(intersected.bed = bed, command = command))
    } else if (return.command == F & return.bed == T) {
        bed = as.data.frame(data.table::fread(outputFileName))
        if (delete.output == T) {unlink(outputFileName)}
        return(bed)
    }

} # END function
