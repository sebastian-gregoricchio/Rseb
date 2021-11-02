#' @title Score matrix NGS data builder at specific regions (by \code{deeptools/computeMatrix} function).
#'
#'
#'
#' @description  This function runs a command line that uses \code{deeptools} to calculate scores per genome regions and to prepare an intermediate file that can be used with \link{plot.density.profile} and \link{plot.density.summary}. Typically, the genome regions are genes, but any other regions defined in a BED file can be used. computeMatrix accepts multiple score files (bigWig format) and multiple regions files (BED format). This tool can also be used to filter and sort regions according to their score.
#'
#'
#'
#' @param mode The type of matrix computation. Allowed values are "reference-point" or "scale-region". No default. \itemize{
#'   \item \code{reference-point}: \cr Reference-point refers to a position within a BED region (e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be plotted;
#'   \item \code{scale-region}: \cr In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.
#'  }
#' @param scoreFileName String vector with the full paths to bigWig file(s) containing the scores to be plotted.
#' @param regionsFileName String vector with the full paths to .BED or .GTF files containing the regions to plot. If multiple bed files are given, each one is considered a group that can be plotted separately. Also, adding a "#" symbol in the bed file causes all the regions until the previous "#" to be considered one group.
#'
#' @param outFileName String containing the full file name to save the gzipped matrix file (.gz) needed by \link{plot.density.profile}.
#' @param outFileNameMatrix If this option is given, then the matrix of values underlying the heatmap will be saved using the indicated name, e.g. IndividualValues.tab. This matrix can easily be loaded into R or other programs. By default \code{NULL}.
#' @param outFileSortedRegions File name in which the regions are saved after skiping zeros or min/max threshold values. The order of the regions in the file follows the sorting order selected. This is useful, for example, to generate other heatmaps keeping the sorting of the first heatmap. Example: Heatmap1sortedRegions.bed. By default \code{NULL}.
#'
#' @param referencePoint Possible choices: TSS, TES, center. The reference point for the plotting could be either the region start (TSS), the region end (TES) or the center of the region. Note that regardless of what you specify, plotHeatmap/plotProfile will default to using "TSS" as the label. By default \code{TSS}.
#' @param nanAfterEnd Logic value. If set (\code{TRUE}), any values after the region end are discarded. This is useful to visualize the region end when not using the scale-regions mode and when the reference-point is set to the TSS. By default \code{FALSE}.
#'
#' @param regionBodyLength Distance in bases to which all regions will be fit. (Default: 1000).
#' @param startLabel Label shown in the plot for the start of the region. Default is TSS (transcription start site), but could be changed to anything, e.g. "peak start". Note that this is only useful if you plan to plot the results yourself and not, for example, with plotHeatmap, which will override this. (Default: "TSS").
#' @param endLabel Label shown in the plot for the region end. Default is TES (transcription end site). See the --startLabel option for more information. (Default: "TES").
#' @param unscaled5prime Number of bases at the 5-prime end of the region to exclude from scaling. By default, each region is scaled to a given length (see the --regionBodyLength option). In some cases it is useful to look at unscaled signals around region boundaries, so this setting specifies the number of unscaled bases on the 5-prime end of each boundary. (Default: 0).
#' @param unscaled3prime Number of bases at the 3-prime end of the region to exclude from scaling. By default, each region is scaled to a given length (see the --regionBodyLength option). In some cases it is useful to look at unscaled signals around region boundaries, so this setting specifies the number of unscaled bases on the 3-prime end of each boundary. (Default: 0).
#'
#' @param upstream Distance upstream of the reference-point selected. (Default: 500).
#' @param downstream Distance downstream of the reference-point selected. (Default: 500).
#' @param binSize Length, in bases, of the non-overlapping bins for averaging the score over the regions length. (Default: 10).
#' @param sortRegions Possible choices: "descend", "ascend", "no", "keep". Whether the output file should present the regions sorted. The default is to not sort the regions. Note that this is only useful if you plan to plot the results yourself and not, for example, with plotHeatmap, which will override this. Note also that unsorted output will be in whatever order the regions happen to be processed in and not match the order in the input files. If you require the output order to match that of the input regions, then either specify "keep" or use computeMatrixOperations to resort the results file. (Default: "keep").
#' @param sortUsing Possible choices: "mean", "median", "max", "min", "sum", "region_length". Indicate which method should be used for sorting. The value is computed for each row.Note that the region_length option will lead to a dotted line within the heatmap that indicates the end of the regions. (Default: "mean").
#' @param sortUsingSamples List of sample numbers (order as in matrix), that are used for sorting by --sortUsing, no value uses all samples, example: --sortUsingSamples 1 3. By default \code{NULL}.
#' @param averageTypeBins Possible choices: "mean", "median", "min", "max", "std", "sum". Define the type of statistic that should be used over the bin size range. (Default: "mean").
#' @param missingDataAsZero Logic value to define if set, missing data (NAs) will be treated as zeros. The default is to ignore such cases (\code{NULL}). If not included, this parameter can be changed later in the function \link{plot.density.profile}.
#' @param skipZeros Logic value to understand whether regions with only scores of zero should be included or not. Default is to include them (\code{FALSE}).
#' @param minThreshold Numeric value. Any region containing a value that is less than or equal to this will be skipped. This is useful to skip, for example, genes where the read count is zero for any of the bins. This could be the result of unmappable areas and can bias the overall results. (Default: \code{NULL}).
#' @param maxThreshold Numeric value. Any region containing a value greater than or equal to this will be skipped. The maxThreshold is useful to skip those few regions with very high read counts (e.g. micro satellites) that may bias the average values. (Default: \code{NULL}).
#' @param blackListFileName A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. (Default: \code{NULL}).
#' @param samplesLabel Labels for the samples. This will then be passed to \link{plot.density.profile} function. The default is to use the file name of the sample. The sample labels should be separated by spaces and quoted if a label itself contains a space E.g. --samplesLabel label-1 "label 2".
#' @param smartLabels Instead of manually specifying labels for the input bigWig and BED/GTF files, this causes deepTools to use the file name after removing the path and extension. (Default: \code{TRUE}).
#' @param scale If set, all values are multiplied by this number. (Default: 1).
#' @param numberOfProcessors Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors. (Default: "max").
#' @param metagene When either a BED12 or GTF file are used to provide regions, perform the computation on the merged exons, rather than using the genomic interval defined by the 5-prime and 3-prime most transcript bound (i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 file is used as input, then columns 2 and 3 are used as an exon. (Default: \code{FALSE}).
#' @param transcriptID When a GTF file is used to provide regions, only entries with this value as their feature (column 3) will be processed as transcripts. (Default: "transcript").
#' @param exonID When a GTF file is used to provide regions, only entries with this value as their feature (column 3) will be processed as exons. CDS would be another common value for this. (Default: "exon").
#' @param transcript_id_designator Each region has an ID (e.g., ACTB) assigned to it, which for BED files is either column 4 (if it exists) or the interval bounds. For GTF files this is instead stored in the last column as a key:value pair (e.g., as ‘transcript_id "ACTB"’, for a key of transcript_id and a value of ACTB). In some cases it can be convenient to use a different identifier. To do so, set this to the desired key. (Default: "transcript_id").
#'
#' @param srun Logic value to define whether the command should be run in \code{srun} mode. By default \code{FALSE}.
#' @param computeMatrix.deeptools.command String to define the command to use to recall the computeMatrix function of deeptools. An example: "/home/user/anaconda3/bin/computeMatrix". By default \code{"/home/USERNAME/anaconda3/bin/computeMatrix"}.
#'
#' @param return.command  Logic value to define whether to return the string corresponding to the command for deeptools. By default \code{FALSE}.
#' @param run.command Logic value to define whether to run the the command line on system terminal and generate the score matrix by deeptools. By default \code{TRUE}.
#'
#' @param quiet Logic value to define if to remove any warning or processing messages. By default \code{FALSE}.
#' @param verbose Logic value to define if to be VERY verbose in the status messages. --quiet will disable this. By default \code{FALSE}.
#'
#'
#'
#' @return The function generates the files indicated by the output parameters. The matrix.gz output file can be read by the function \link{read.computeMatrix.file}.
#'
#'
#'
#' @examples
#' computeMatrix.deeptools(
#'    mode = "reference-point",
#'    scoreFileName = c("path_to/signal_file1.bw", "path_to/signal_file2.bw"),
#'    regionsFileName = c("path.to/regions1.bed", "path.to/regions2.bed"),
#'    upstream = 1000,
#'    downstream = 1000,
#'    outFileName = "path_to/output_matrix.gz",
#'    computeMatrix.deeptools.command = "/home/user/anaconda3/bin/computeMatrix",
#'    referencePoint = "peakMax")
#'
#' computeMatrix.deeptools(
#'    mode = "scale-regions",
#'    scoreFileName = c("path_to/signal_file1.bw", "path_to/signal_file2.bw"),
#'    regionsFileName = c("path.to/regions1.bed", "path.to/regions2.bed"),
#'    upstream = 1000,
#'    downstream = 1000,
#'    regionBodyLength = 300,
#'    startLabel = "geneStart",
#'    endLabel = "geneEnd",
#'    outFileName = "path_to/output_matrix.gz",
#'    computeMatrix.deeptools.command = "/home/user/anaconda3/bin/computeMatrix",
#'    referencePoint = "peakMax")
#'
#'
#'
#' @details To know more about the deeptools's \code{computeMatrix} function see the package manual at the following link: \cr \url{https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html}.
#'
#'
#'
#' @export computeMatrix.deeptools



computeMatrix.deeptools =
  function(
    # basic input parameters
    mode,
    scoreFileName,
    regionsFileName,

    # output parameters
    outFileName,
    outFileNameMatrix = NULL,
    outFileSortedRegions = NULL,

    # seference-point mode
    referencePoint = "TSS",
    nanAfterEnd = FALSE,

    # scale-region mode
    regionBodyLength = 1000,
    startLabel = "TSS",
    endLabel = "TES",
    unscaled5prime = 0,
    unscaled3prime = 0,

    # Common parameters
    upstream = 500,
    downstream = 500,
    binSize = 10,
    sortRegions = "keep",
    sortUsing = "mean",
    sortUsingSamples = NULL,
    averageTypeBins = "mean",
    missingDataAsZero = FALSE,
    skipZeros = FALSE,
    minThreshold = NULL,
    maxThreshold = NULL,
    blackListFileName = NULL,
    samplesLabel = NULL,
    smartLabels = TRUE,
    scale = 1,
    numberOfProcessors = "max",
    metagene = FALSE,
    transcriptID = "transcript",
    exonID = "exon",
    transcript_id_designator = "transcript_id",

    # path to deeptools
    srun = FALSE,
    computeMatrix.deeptools.command = paste0("/home/", Sys.getenv("USERNAME"), "/anaconda3/bin/computeMatrix"),

    # others parameters for this function in R
    return.command = FALSE,
    run.command = TRUE,

    # Verbose options
    quiet = FALSE,
    verbose = FALSE
  ) {

    #-----------------------------#
    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)   #
    #-----------------------------#

    ######################################################################################
    # Create function to add single quote to string variables
    add.quotes = function(x) {return(sapply(x, function(x){paste("'", x, "'", sep = "")}, USE.NAMES = F))}



    ###### Generation of the string to run the command on bash
    if (srun == TRUE) {
      command = paste("srun", computeMatrix.deeptools.command)
    } else {command = computeMatrix.deeptools.command}



    ###### Check and save the mode of computation
    if (!(mode %in% c("reference-point", "scale-regions"))) {
      return(warning("The mode for the computation of the matrix must be one among: 'reference-point', 'scale-regions'."))
    } else {command = paste(command, mode)}



    ###### Add the bw files and the regions bed/gtf
    command = paste(command,
                    "--scoreFileName", paste(add.quotes(scoreFileName), collapse = " "),
                    "--regionsFileName", paste(add.quotes(regionsFileName), collapse = " "))



    ###### Add output arguments
    command = paste(command, "--outFileName", outFileName)
    if (!is.null(outFileNameMatrix)) {command = paste(command, "--outFileNameMatrix", add.quotes(outFileNameMatrix))}
    if (!is.null(outFileSortedRegions)) {command = paste(command, "--outFileSortedRegions", add.quotes(outFileSortedRegions))}



    ###### Add mode-dependent parameters
    # REFERENCE-POINT parameters
    if (mode == "reference-point") {
      command = paste(command, "--referencePoint", add.quotes(referencePoint))
      if (nanAfterEnd == T) {command = paste(command, "--nanAfterEnd")}

    } else {
      # SCALE-REGIONS parameters
      command =
        paste(
          command,
          "--regionBodyLength", regionBodyLength,
          "--startLabel", add.quotes(startLabel),
          "--endLabel", add.quotes(endLabel),
          "--unscaled5prime", unscaled5prime,
          "--unscaled3prime", unscaled3prime)
    }



    ###### Add common fixed-values
    command =
      paste(command,
            "--upstream", upstream,
            "--downstream", downstream,
            "--binSize", binSize,
            "--sortRegions", sortRegions,
            "--sortUsing", sortUsing,
            "--averageTypeBins", averageTypeBins,
            "--numberOfProcessors", numberOfProcessors,
            "--transcriptID", transcriptID,
            "--exonID", exonID,
            "--transcript_id_designator", transcript_id_designator,
            "--scale", scale)



    ###### Add common optional parameters
    if (!is.null(sortUsingSamples)) {command = paste(command, "--sortUsingSamples", paste(sortUsingSamples, collapse = " "))}
    if (!is.null(minThreshold)) {command = paste(command, "--minThreshold", paste(minThreshold, collapse = " "))}
    if (!is.null(maxThreshold)) {command = paste(command, "--maxThreshold", paste(maxThreshold, collapse = " "))}
    if (!is.null(blackListFileName)) {command = paste(command, "--blackListFileName", add.quotes(blackListFileName))}
    if (!is.null(samplesLabel)) {command = paste(command, "--samplesLabel", paste(add.quotes(samplesLabel), collapse = " "))}

    if (missingDataAsZero == T) {command = paste(command, "--missingDataAsZero")}
    if (skipZeros == T) {command = paste(command, "--skipZeros")}
    if (smartLabels == T) {command = paste(command, "--smartLabels")}
    if (metagene == T) {command = paste(command, "--metagene True")}

    if (quiet == T) {command = paste(command, "--quiet")}
    if (verbose == T) {command = paste(command, "--verbose")}



    ###### Run and/or export the command
    if (run.command == T) {
      system(command)
      message("Generation of the output file(s): done.")}

    if (return.command == T) {return(command)}

} # END function
