#' @title Script generator for Integrative Genomics Viewer (IGV) batch tasks.
#'
#' @description Helps to the generation of a script file that can be run on IGV to generate multiple screenshots at specific genomic regions.
#'
#' @param loci_vector Either a gene name vector (e.g. \code{c("Gapdh", "Spi1", ...)}) or a regions vector (eg. \code{c('chr1:253000-256503', ...)}. All IGV formats are allowed.
#' @param input_type Define the input type. Allowed values are \code{genes} and \code{regions}.
#' @param biomart Defines the \code{biomart} parameter for \code{biomaRt} package, by default \code{ensembl}.
#' @param dataset Defines the \code{dataset} parameter for \code{biomaRt} package, by default \code{mmusculus_gene_ensembl}.
#' @param reference_genome [optional] Defines the genome to use, e.g. "mm10", "hg19", ... . By default \code{NULL}.
#' @param fivePrime Numeric value to define how many bases [bp] exapand from full gene position at it's 5'-end, default 1000bp.
#' @param threePrime Numeric value to define how many bases [bp] exapand from full gene position at it's 3'-end, default 1000bp.
#' @param snap_names [optional] String vector to define the names of images (without extention), by default uses \code{loci_vector}.
#' @param IGV_batch_file String for the batch_script_file_name/path, by default \code{<working_directory>/IGV_batch.txt}.
#' @param snap_image_format String to define the format of the images, e.g. "png", "jpeg", "svg", ... . By default \code{png}.
#' @param snap_directory String for the output directory for the snapshoots. By default <working_directory>.
#' @param maxPanelHeight Numeric value to define the height in pixel of the IGV pannel that will be captured on IGV.
#' @param session [optional] FULL path to an IGV session file (session.xml) to use for the images. By default \code{NULL}.
#' @param exit Logical value to indicate whether exit IGV after image capture ended. By default \code{FALSE}.
#' @param help Logical value to indicate whether display the help. By default \code{FALSE}.
#'
#' @return Exports a .txt file ready-to-use on IGV.
#'
#' @details For more info on how batch tasks work on IGV see: \cr \url{https://software.broadinstitute.org/software/igv/PortCommands}.
#'
#' @export IGVsnap
#'
#' @import biomaRt
#' @import dplyr

######################
## IGVsnap function ##
######################
IGVsnap = function(loci_vector,
                   input_type, # 'genes' or 'regions'
                   biomart = "ensembl",
                   dataset = "mmusculus_gene_ensembl",
                   reference_genome = NULL,
                   fivePrime = 1000, #bp
                   threePrime = 1000,
                   snap_names = NULL,
                   IGV_batch_file = paste(getwd(), "/IGV_batch.txt", sep = ""),
                   snap_image_format = "png",
                   snap_directory = getwd(),
                   maxPanelHeight = 1000,
                   session = NULL,
                   exit = FALSE,
                   help = FALSE) {

  # Install packages from bioconductor
  pkg = "biomaRt"
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    if(!require(pkg, character.only = TRUE)) stop(paste(pkg,"package not found."))}

  # check parameters
    help_message = c(
    "Help for 'IGVsnap' function from 'Rseb' package: \n", "\n",
    "      loci_vector   either a gene name vector [e.g. c('Gapdh', 'Spi1', ...)] or a regions vector [eg. c('chr1:253000-256503', ...). All IGV formats are allowed. \n",
    "       input_type   define whether 'genes' or 'regions'. \n",
    "          biomart   'biomart' parameter for biomaRt package, by default 'ensembl'. \n",
    "          dataset   'dataset' parameter for biomaRt package, by default 'mmusculus_gene_ensembl'. \n",
    " reference_genome   [optional] to define the genome to use 'mm10', 'hg19, ..., 'NULL' by default. \n",
    "        fivePrime   how many bases [bp] exapand from full gene position at it's 5'-end, default 1000bp. \n",
    "       threePrime   how many bases [bp] exapand from full gene position at it's 3'-end, default 1000bp. \n",
    "       snap_names   [optional] names of images (without extention), by default uses 'loci_vector'. \n",
    "   IGV_batch_file   batch_script_file name/position, by default <working_directory>/IGV_batch.txt. \n",
    "snap_image_format   image format such as 'png', 'jpeg', 'svg', ..., by default 'png'. \n",
    "   snap_directory   output directory for snapshoots, by default <working_directory>. \n",
    "   maxPanelHeight   heigth in pixel of the IGV pannel that will be captured. \n",
    "          session   [optional] define an IGV session file .xml to use for the images, 'NULL' by default - USE FULL PATH -. \n",
    "             exit   logical to indicate whether exit IGV after image capture ended, FALSE by default. \n",
    "             help   logical to indicate whether display the help, FALSE by default. \n",
    "\n", "More info at https://software.broadinstitute.org/software/igv/PortCommands")

  if (exists("input_type") == F | exists("loci_vector") == F | help == T | (!(input_type %in% c("genes", "regions")))) {return(message(help_message))}

  # Retrieves regions
  if (input_type == "genes") {
    loci_vector = sort(unique(loci_vector))

    require(biomaRt)
    # to list all the datasets availables
    # listDatasets(useMart("ensembl"))

    # for getBM():
    # EXAMPLES
    # mm10 = useMart("ensembl", dataset="mmusculus_gene_ensembl")
    # listFilters(mm10)     /to filter output
    # listAttributes(mm10)  /output columns

    genome = useMart(biomart = biomart,
                     dataset = dataset)

    # loading gene list
    require(dplyr)
    gene_positions = getBM(mart = genome, # mart object, see "useMart" function
                           attributes = c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name"), # columns in output
                           filters = "external_gene_name", # select only certain values,
                           values = loci_vector) %>% # values for the filter
    arrange(factor(external_gene_name, levels = loci_vector))

    gene_positions$start_position = gene_positions$start_position-fivePrime
    gene_positions$end_position = gene_positions$end_position+threePrime

    IGV_positions = data.frame(position = paste("chr", gene_positions$chromosome_name, ":",
                                                gene_positions$start_position, "-",
                                                gene_positions$end_position,
                                                sep = ""),
                               gene = gene_positions$external_gene_name)

    if (is.null(snap_names)) {snap_names = IGV_positions$gene}
    list = IGV_positions$position

    } else if (input_type == "regions") {
      if (is.null(snap_names)) {snap_names = loci_vector}
      list = unique(loci_vector)
    }

  # Initialization/Creation of the batch_file
  write(file = IGV_batch_file,
        x = paste("# IGV_batch_file,", as.character(date())))

  # Insert loading exiting session if required
  if (!is.null(session)) {
    write(file = IGV_batch_file,
          x = paste("load", session),
          append = T)
  } else if (!is.null(reference_genome)){
    write(file = IGV_batch_file,
          x = paste("genome", reference_genome),
          append = T)
  }


  # Set pannel size (height in pixel)
  write(file = IGV_batch_file,
        x = paste("maxPanelHeight", maxPanelHeight),
        append = T)


  # Insert output directory for the snapshots
  write(file = IGV_batch_file,
        x = paste("snapshotDirectory", snap_directory),
        append = T)

  for (i in 1:length(list)) {
    # Write the position where to go
    write(file = IGV_batch_file,
          x = paste("goto", list[i]),
          append = T)

    # Write the command to take a snapshot with the name of the gene/position with the chosen extension
    write(file = IGV_batch_file,
          x = paste("snapshot ", snap_names[i], "_snapshot.", snap_image_format, sep =""),
          append = T)
  }

  # Add command to go to whole genome (to indicate that finished)
  write(file = IGV_batch_file,
        x = "goto all",
        append = T)

  # Exit command if required
  if (exit == TRUE) {
    write(file = IGV_batch_file,
          x = "exit",
          append = T)
  }

  message("The following batch script have been generated: \n",
          IGV_batch_file,
          "\n", "\n",
          "The final snapshoot will be generated in the following folder: \n",
          snap_directory)

} # END funtion
