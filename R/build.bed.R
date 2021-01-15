#' @title Bed generator
#'
#' @description Function that helps the building of a bed file providing the columns. It enables also the specification of the track line for software such as IGV in order to pre-define colors, track name, etc.
#'
#' @param chr String vector containing the name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
#' @param start Numeric vector indicating the starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#' @param end Numeric vector indicating the ending position of the feature in the chromosome or scaffold.
#' @param name String vector defining the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode. If set as \code{NULL} (default) and the column is required, the names will correspond to the mid-point of the region.
#' @param score A single value or a numeric vector with a score between 0 and 1000. If the track line \code{useScore} attribute is set as \code{TRUE} for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). By default 0.
#' @param strand A single character or a string vector defining the strand: either "." (=no strand) or "+" or "-". By default ".".
#' @param thickStart A numeric vector indicating the starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part (default value, \code{thickStart = NULL}) it will be used the \code{start} value.
#' @param thickEnd A numeric vector indicating the ending position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part (default value, \code{thickStart = NULL}) it will be used the \code{end} value.
#' @param itemRgb A single value or a string vector containing the colors for each feature. It can be expressed as an RGB value of the form R,G,B (e.g. "255,0,0") or as any other R-supported color name (it will be converted automatically to RGB version). By default \code{NULL}. If the track line \code{itemRgb.ON} attribute is set as \code{TRUE}, this color value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
#' @param blockCount A single number or a numeric vector indicating the number of blocks (exons) in the BED line. By default \code{NULL}.
#' @param blockSizes A vector containing a comma-separated list of the block sizes. The number of items in this list should correspond to \code{blockCount}. By default \code{NULL}.
#' @param blockStarts A vector containing a comma-separated list of block starts. All of the \code{blockStart} positions should be calculated relative to \code{start}. The number of items in this list should correspond to \code{blockCount}. By default \code{NULL}.
#'
#' @param track.name A string defining the track label that will be displayed to the left of the track in the Genome Browser window, and also the label of the track control at the bottom of the screen. The name can consist of up to 15 characters. It is recommended that the track_label be restricted to alpha-numeric characters and spaces to avoid potential parsing problems. By default \code{NULL}.
#' @param display.mode A string that defines the initial display mode of the annotation track. Values for \code{display.mode} include: "hide", "dense", "full", "pack", "squish". By default \code{NULL}.
#' @param itemRgb.ON Logic value to define whether this attribute should be set to "On", the Genome Browser will use the RGB value shown in the \code{itemRgb} field in each data line of the associated BED track to determine the display color of the data on that line. If the \code{itemRgb} values are not provided, this parameter will be ignored. By default \code{TRUE}.
#' @param useScore Logic value to define if the \code{score} field in each of the track's data lines should be used to determine the level of shading in which the data is displayed. By default \code{FALSE}.
#' @param colorByStrand A vector composed by two strings for two colors, either in RGB comma separated format (eg. "0,250,30") or any R-supported color string (they will be converted automatically to RGB format). The order of color sets is c("strand +", "strand -"). Parameter ignored when \code{itemRgb} is active/provided. By default \code{NULL}.
#' @param track.base.color A single string defining the main color for the annotation track. The track color consists of three comma-separated RGB values from 0-255 (eg. "0,250,30") or any R-supported color string (it will be converted automatically to RGB format). Parameter ignored when \code{itemRgb} or \code{colorByStrand} are active/provided. By default \code{NULL}.
#'
#' @param sort Logic value to define whether to sort the bed using the function \link{sort.bed}. By default \code{TRUE}.
#' @param bed.file.name If a string with a full path to a bed_file is provided, the function will export the bed as a txt file. By default \code{NULL}.
#' @param export.track.line Logic value to define if the track line should be exported. When \code{bed.file.name = NULL} this parameter is ignored. By default \code{TRUE}.
#' @param return.data.frame Logic value to define if the to return the data.frame corresponding to the bed (it will show the columns names). By default \code{FALSE}.
#' @param force.generation Force the generation of bed even when certain errors occur (eg. score > 1000, start > end). By default \code{FALSE}.
#'
#'
#' @return If required the function can export a bed file with or without the track line, return a data.frame (with column names) corresponding to the bed generated, or both. The bed file could be automatically sorted settin the parameter \code{sort = TRUE}.
#'
#' @references
#' \itemize{
#'    \item More information about bed format are available at the following link: \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}.
#'    \item More information about track line parameters are available at the following link: \url{https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#lines}.
#' }
#'
#' @export build.bed


build.bed = function(
  # bed columns
  chr,
  start,
  end,
  name = NULL,
  score = 0,
  strand = ".",
  thickStart = NULL,
  thickEnd = NULL,
  itemRgb = NULL,
  blockCount = NULL,
  blockSizes = NULL,
  blockStarts = NULL,

  # track lines
  track.name = NULL,
  display.mode = NULL,
  itemRgb.ON = T,
  useScore = F,
  colorByStrand = NULL,
  track.base.color = NULL,

  # export parameters
  sort = T,
  bed.file.name = NULL,
  export.track.line = TRUE,
  return.data.frame = F,
  force.generation = F
  ) { # BEGIN function

  # -------------------------------------------------------------------------- #

  # Check size of vectors for bed generation
  ##############################################################################
  row.number = length(chr)

  param.check = function(param, row.number){
    ifelse(test = if (!is.null(param)) {length(param) > 1},
           yes = length(param) == row.number,
           no = TRUE)
  }

  length.check = c(
    length(start) == row.number,
    length(end) == row.number,
    param.check(name, row.number),
    param.check(score, row.number),
    param.check(strand, row.number),
    param.check(thickStart, row.number),
    param.check(thickEnd, row.number),
    param.check(itemRgb, row.number),
    param.check(blockCount, row.number),
    param.check(blockSizes, row.number),
    param.check(blockStarts, row.number)
  )

  if (length(unique(length.check)) > 1) {return(warning(paste("The expected number of elements in each colmun/parameter must be", row.number, "or 1.")))}


  # Check that START is before END and force if required
  if (-1 %in% sign(end - start) & force.generation == F) {
   return(warning("At least one value in END vector is lower than the relative START value. If you want to force the bed generation set the parameter 'force.generation' as TRUE."))
  }

  # Check that that scores are between 0 and 1000
  if (-1 %in% sign(1000 - score) & force.generation == F) {
    return(warning("At least one value in SCORE vector is grater than 1000 (allowed values between 0-1000). If you want to force the bed generation set the parameter 'force.generation' as TRUE."))
  }

  # Convert non-RGB colors (if required)
  if (!is.null(itemRgb)) {
    if (!grepl(",", itemRgb[1])) {
      itemRgb = sapply(itemRgb, function(x){paste(as.vector(col2rgb(x)), collapse = ",")}, USE.NAMES = F)
    }
  }


  # Names generation
  if (is.null(name)) {name = round((start + end)/2)}

  # Define thickStart and thickEnd if RGB is defined
  if (!is.null(itemRgb)) {
    if (is.null(thickStart)) {thickStart = start}
    if (is.null(thickEnd)) {thickEnd = end}
  }

  # Check blocks
  blocks.check = c(!is.null(blockCount), !is.null(blockSizes), !is.null(blockStarts))
  if (length(unique(blocks.check)) > 1) {
    return(warning("If you want to define blocks, you need to provide all the 3 vectors for 'blockCount', 'blockSizes', 'blockStarts'."))
  }

  # Check track label length
  if (!is.null(track.name)) {
    if (nchar(track.name) > 15 | length(track.name) > 1) {return(warning("The label of the track.name must be a string with a number of characters <=15."))}
  }


  # TRACK LINE BUILDING
  ##############################################################################
  track.line = NULL
  if (!is.null(bed.file.name)) {
    if (export.track.line == TRUE) {
      if (!is.null(track.name)) {if (length(track.name) == 1) {track.line = c(track.line, paste('track name="', track.name, '"', sep = ""))}}
      if (!is.null(display.mode)) {if (length(display.mode) == 1 & display.mode %in% c("hide", "dense", "full", "pack", "squish")) {track.line = c(track.line, paste('visibility="', display.mode, '"', sep = ""))}}
      if (itemRgb.ON == TRUE & !is.null(itemRgb)) {track.line = c(track.line, 'itemRgb="on"')}
      if (useScore == TRUE & length(score) > 1) {track.line = c(track.line, 'useScore="1"')}

      if (itemRgb.ON == FALSE) {
        if (!is.null(colorByStrand)) {
          if (length(colorByStrand) == 2) {
            if (!grepl(",", colorByStrand[1])) {strand.colors = paste(sapply(colorByStrand, function(x){paste(as.vector(col2rgb(x)), collapse = ",")}, USE.NAMES = F), collapse = " ")}
            track.line = c(track.line, paste('colorByStrand="', strand.colors, '"', sep = ""))
            }
          } # END if for colorByStrand

        if (is.null(colorByStrand) & !is.null(track.base.color)) {
          if (length(track.base.color) == 1) {
            if (!grepl(",", track.base.color)) {base.color = paste(as.vector(col2rgb(track.base.color)), collapse = ",")}
            track.line = c(track.line, paste('color="', base.color, '"', sep = ""))
          } # END if for track.base.color
        }
      } # END if for itemRGB OFF

      track.line = paste(track.line, collapse = " ")
    } # END if for export line
  } # END if for bed_file

  # DATA.FRAME BUILDING
  ##############################################################################
  # Select the columns required for the bed
  columns_list = list(chr = chr,
                      start = start,
                      end = end,
                      name = name,
                      score = score,
                      strand = strand,
                      thickStart = thickStart,
                      thickEnd = thickEnd,
                      itemRgb = itemRgb,
                      blockCount = blockCount,
                      blockSizes = blockSizes,
                      blockStarts = blockStarts)

  columns.to.keep = c()
  for (i in 1:length(columns_list)) {
    if (!is.null(columns_list[[i]])) {columns.to.keep = c(columns.to.keep, i)}
  }

  # build data.frame and sort if required
  bed = data.frame(columns_list[columns.to.keep])
  if (sort == T) {bed = Rseb::sort.bed(bed)}


  # OUTPUT GENERATION
  ##############################################################################
  # Coerce track.line when empty
  if (is.null(track.line)) {track.line = ""}

  # Export file if required
  if (!is.null(bed.file.name)) {
    if (length(bed.file.name) == 1) {
      if (export.track.line == T & track.line != "") {
          write(x = track.line, file = bed.file.name)
          write.table(x = bed, file = bed.file.name,
                      quote = F, sep = "\t", row.names = F, col.names = F, append = T)
          message(paste("Bed file exported as -> ", bed.file.name, sep = ""))
      } else {
        write.table(x = bed, file = bed.file.name,
                    quote = F, sep = "\t", row.names = F, col.names = F)
        message(paste("Bed file exported as -> ", bed.file.name, sep = ""))
      }
    }
  }

  # Return data.frame if required
  if (return.data.frame == T) {return(bed)}

} # END function
