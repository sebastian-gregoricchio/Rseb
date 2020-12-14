#' @title Generator of a bed file for enzymatic restriction sites.
#'
#' @description The function allows to create a bed file that can be added on IGV both as regions and track. It will show the restriction sites of a sequences if starting from the cut positions depending on sequence lenght. Chromosome, start and end of the input sequence are required.
#'
#' @param cut_positions A numeric vector with the list of the restriction/cut positions.
#' @param chromosome Chromosome number of the region analyzed.
#' @param genome_start Start position on the genome of the region analyzed.
#' @param return_bed Logic value to define if to return the bed as data.frame. By default \code{TRUE}.
#' @param export_bed_file Logic value to define if to export the resulting .bed file. By default \code{FALSE}.
#' @param output_file_name String corresponding to the path to the exported .bed file. By default \code{"<working.directory>/restriction_positions.bed"}.
#' @param enzyme_cut_length Numeric value to define the length of cut of the restriction enzyme. By default 4.
#' @param include_region_description Logic value to define whether to include a fourth column containing the region name define by the parameter \code{region_description}. By default \code{TRUE}.
#' @param region_name Regions base name. Automatically it will be added a number to the base name. By default \code{"site"}, the resulting regions will be: site_1, site_2, ... .
#' @param append Logic value to define if to append the result to the file. By default \code{FALSE}, the file will be overwritten.
#'
#' @return If required, it will be returned a classic bed file (chr start end [name]) with the regions centered on the cut position in the genome.
#'
#' @examples
#' restriction.sites.to.bed(cut_positions = c(230, 235, 1250, 36),
#'                          chromosome = 10,
#'                          genome_start = 1205126,
#'                          region_name = "EcoRI_cut_site")
#'
#' @details
#' To map the positions of restriction enzymes it is possible to use \url{http://restrictionmapper.org/} with the option \code{Map (version 3)}.
#'
#' @export restriction.sites.to.bed

restriction.sites.to.bed = function(
  cut_positions, #must be a vector with the list of the positions eg. c(230, 235, 1250, 36)
  chromosome,    #number of the chromosome
  genome_start,  #start in bp
  return_bed = TRUE,
  export_bed_file = FALSE,
  output_file_name = paste(getwd(), "restriction_positions.bed", sep = "/"),
  enzyme_cut_length = 4,
  include_region_description = TRUE,
  region_name = "site",
  append = FALSE)

{

# Sorting the positions vector
cut_positions = sort(cut_positions, decreasing = FALSE)

#Creation of the variables
chr = c()
start = c()
end = c()
regions = c()

# It assigns 1 to enzyme-cut positions and 0 where there is not
for (i in c(1:length(cut_positions))) {
  chr[i] = paste("chr", chromosome, sep="", collapse="")
  start[i] = cut_positions[i] + genome_start
  end[i] = start[i] + enzyme_cut_length
  regions[i] = paste(region_name, (i-1), sep="_", collapse = "")}


# Removing of all positions that do not contain a cut site, or 'keep all values different from 0' for each vector
# chr = chr[chr != 0]
# start = start[start != 0]
# end = end[end != 0]

# Creation of the table/bed combining the vectors in the right order
if (include_region_description == TRUE)
(bed_file_cut_positions = data.frame(chr, start, end, regions)) else (
  bed_file_cut_positions = data.frame(chr, start, end))

# Exporting the .bed file
if (export_bed_file == TRUE) {
  write.table(x = bed_file_cut_positions,
              file = output_file_name,
              sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE,
              append = append)

  message("The .bed file has been saved as ", output_file_name)
}

#check the results
if (return_bed == TRUE) {return(bed_file_cut_positions)}
}
