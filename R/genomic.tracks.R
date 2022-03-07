#' @title Genomic tracks plotter
#'
#' @description The functions allows to plot different types of genomic data (bigWig, bed, bedpe) at a specific genomic region. It is possible to highlight specific regions and the gene annotations are plotted automatically at the bottom of all the tracks.
#'
#' @param tracks A vector indicating the list of full paths of the files/tracks/signals to plot. Supported formats: bed/bd/narrowPeak/broadPeak, bw/bigWig/bigwig, bedpe.
#' @param genomic.region An atomic string indicating the genomic region into which restrict the final plot in the format 'chr1:1234-5678'.
#' @param genome An atomic string indicating the genome to use for the annotations. Allowed values are:
#' \itemize{
#'   \item \code{hg19}: loads an 'EnsDb' object from the library \code{EnsDb.Hsapiens.v75};
#'   \item \code{hg38}: loads an 'EnsDb' object from the library \code{EnsDb.Hsapiens.v86};
#'   \item \code{mm10}: loads an 'EnsDb' object from the library \code{EnsDb.Mmusculus.v79};
#'   \item \emph{custom 'EnsDb' object}: provide an 'EnsDb' object manually generated; visit the page \url{https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#102_building_annotation_packages} for more information.
#'  }
#' @param track.labels A vector indicating the labels to use for each track (genome annotation track excluded). By default \code{NULL}: the file base-name will be used.
#' @param track.labels.fontzise A numerical value to indicate the font size of the track labels. Default value \code{5}.
#' @param track.labels.position A two-element numeric vector passed to \code{xlim} function for the the definition of the frame size of the track labels. Default value \code{c(-0.1, 0)}.
#' @param track.colors A string vector indicating the color to use for each track (genome annotation track excluded). If only one value is provided it will be used for all the tracks. Default value \code{"#000000"} ("black").
#' @param grouping A single numerical vector or a list of numeric vectors. Each list's element indicates the indexes corresponding to the tracks (1 = first track, 2 = second track, etc) for which the y-axes should be normalized. Each element will be taken into account in the order. Default value \code{NULL}.
#' @param gene.annotation.color A string indicating the color to use for the genome annotation track.
#' @param expand.bed A logical value to define whether overlapping regions in a bed should be plotted on different levels. Default \code{TRUE}.
#' @param arcs.direction A string indicating the direction on which arcs should be plotted for bedpe files. Available options \code{"up"} or \code{"down"}. Default value \code{"down"}.
#' @param fraction.arc.base A numerical value indicating the fraction of total plot height to be used as arc base thickness. By default \code{0.025} (2.5\% of the track height).
#' @param highlight.bed Either a string indicating the full path to a bed file or a data.frame in BED3 format (chr, start, end) containing regions that should be highlighted in the plot. Regions included in the genomic range will be automatically selected. By default \code{NULL}.
#' @param highlight.color A string indicating the color to use for the regions to highlight in the plot. By default \code{'yellow'}.
#' @param highlight.transparency A numerical value indicating the transparency (alpha) to use for the highlighted regions. Default value \code{0.15}.
#' @param missing.data.as.zero.bw A logical value to define wthere missing data in the bigWigs should be converted to zeros. Default \code{FALSE}.
#' @param smooth.bigWig.signal Logical value to indicate whether the bigWig signals should be smoothed (by loess x ~ y function). By default \code{TRUE}.
#' @param smooth.bigWig.loess.span Numerical value to indicate the span value for the loess function used to smooth bigWig signals. By default \code{0.05}.
#' @param plot.bigWig.area Logical value to indicate whether the bigWig profile should be filled or not. If \code{FALSE} only the signal outline will be plotted. By default \code{TRUE}.
#' @param bigWig.range.label.size A numerical value to indicate the font size of the bigWig signal range. Default value \code{2.5}.
#' @param score.bed.shadow Logical value to define whether the filling intensity of the bed segments should reflect the score of each signal. By default \code{FALSE}.
#' @param height.ratios Numerical vector of relative track heights, passed to 'rel_heights' parameter of \code{cowplot::plot_grid()}. For example, in a two-row grid, rel_heights = c(2, 1) would make the first column twice as wide as the second column. Value \code{1} indicates that all the tracks should have the same size. By default \code{NULL}, automatic ratios will be computed by this function.
#' @param width.ratios Numerical vector of relative labels vs tracks widths, passed to 'rel_widths' parameter of \code{cowplot::plot_grid()}. For example, in a two-column grid, rel_widths = c(2, 1) would make the first column twice as wide as the second column. Value \code{1} indicates that all the tracks should have the same size. By default \code{c(1,5)} (1 label : 5 tracks).
#'
#' @return The function returns a named list containing:
#' \itemize{
#'   \item \code{configuration}: data.frame with the parameters used to build the plot(s);
#'   \item \code{highlighted.region}: data.frame with the regions used for the highlighting;
#'   \item \code{single.track.list}: a named list containing each single track plot used for the creation of the multi.track.plot;
#'   \item \code{single.label.plot.list}: a named list containing each single track label plot used for the creation of the multi.track.plot;
#'   \item \code{multi.track.plot}: the assembled multi.track labelled plot.
#'  }
#'
#'
#' @export genomic.tracks
#'
#'
#'
genomic.tracks =
  function(
    tracks,
    genomic.region,
    genome,
    track.labels = NULL,
    track.labels.fontzise = 5,
    track.labels.position = c(-0.1, 0),
    track.colors = "#000000",
    grouping = NULL,
    gene.annotation.color = "darkblue",
    expand.bed = TRUE,
    arcs.direction = "down",
    fraction.arc.base = 0.025,
    highlight.bed = NULL,
    highlight.color = "yellow",
    highlight.transparency = 0.15,
    missing.data.as.zero.bw = FALSE,
    smooth.bigWig.signal = TRUE,
    smooth.bigWig.loess.span = 0.05,
    plot.bigWig.area = TRUE,
    bigWig.range.label.size = 2.5,
    score.bed.shadow = FALSE,
    height.ratios = NULL,
    width.ratios = c(1,5)) {

    #-----------------------------#
    # Check if Rseb is up-to-date #
    Rseb::actualize(update = F, verbose = F)   #
    #-----------------------------#


    ################### LIBRARIES ###################
    require(dplyr)
    require(ggbio)
    require(ggforce)
    require(ggplot2)
    require(GenomicRanges)
    #################################################




    # General controls and check of inputs
    ### Check genomic region and generate a bed file
    genomic.region = gsub(pattern = "[,]", replacement = "", x = genomic.region)
    genomic.region = gsub(pattern = "[-]|[+]|[:]", replacement = "_", x = genomic.region)
    genomic.region = strsplit(genomic.region, "_")[[1]]
    if (genomic.region[3] <= genomic.region[2]) {return(warning("The genomic region to be shown must be in the format chr1:12345-67890, where the end position must be greater than the start position and the window must have a length > 1."))}


    ### Genome check
    if (class(genome) == "character") {
      genome = tolower(genome)
      if ((!(genome %in% c("hg19", "hg38", "mm10")))) {return(warning("Genome not recognized. The supported genomes are: 'Hg19', 'Hg38', 'Mm10'. Alternatively a manually build EnsDB file can be provided (see: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#102_building_annotation_packages)."))}
      if (genome == "hg19") {ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75}
      if (genome == "hg38") {ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86}
      if (genome == "mm10") {ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79}
    } else if ("EnsDb" %in% class(genome)) {
      ensdb = genome
    } else {return(warning("Genome not recognized. The supported genomes are: 'Hg19', 'Hg38', 'Mm10'. Alternatively a manually build EnsDB file can be provided (see: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#102_building_annotation_packages)."))}


    ### Auto-identification of the tracks format
    tracks_format = sapply(tracks,
                           function(x){
                             if (grepl("[.]bigwig$|[.]bw$|[.]bigWig$", x, ignore.case = T)) {return("bigWig")
                             } else if (grepl("[.]bed$|[.]bd$|[.]narrowPeak$|[.]broadPeak$", x, ignore.case = T)) {return("bed")
                             } else if (grepl("[.]bedpe$", x, ignore.case = T)) {return("arcs")
                             } else {return(return(warning(paste0("The format of the track file '", x, "' is not recognized. Available formats are: .bw|bigWig|bw, .bed|bd|narrowPeak|broadPeak, .bedpe."))))}
                           },
                           USE.NAMES = F)



    ### Track labels
    if (!is.null(track.labels)) {
      if (length(track.labels)!=1 & ((length(track.labels) < length(tracks)) | (length(track.labels) > length(tracks)))) {
        return(warning("The number of track labels does not correspond to the number of tracks provided."))}
    } else if (is.null(track.labels)) {
      track.labels = gsub(pattern = "[.].*$", replacement = "", x = basename(tracks))
    }



    ### Grouping check
    if (!is.null(grouping) & class(grouping) != "list") {
      grouping = list(grouping)
    }


    ### Track colors
    if (!is.null(track.colors)) {
      if (length(track.colors)!=1 & ((length(track.colors) < length(tracks)) | (length(track.colors) > length(tracks)))) {
        message("The number of track colors does not correspond to the number of tracks provided. The track colors is coerced to 'black' for all tracks.")
        track.colors = "black"}
    } else if (is.null(track.colors)) {
      track.colors = "black" }

    if (F %in% Rseb::is.color(track.colors)) {
      return(warning(paste0("The following track.colors values are not colors: ", paste(track.colors[!Rseb::is.color(track.colors)], collapse = ", "), ".")))}


    ### Arcs direction
    if (tolower(arcs.direction) %in% c("up", "down")) {
      arcs.direction = ifelse(test = tolower(arcs.direction) == "up",
                              yes = 1,
                              no = -1)
    } else {return(warning("The arcs direction must be 'up' or 'down'."))}



    ### Read highlight region
    if (!(is.null(highlight.bed))) {
      if ("data.frame" %in% class(highlight.bed)) {colnames(highlight.bed)[1:3] = c("chr", "start", "end")
      } else {
        highlight.bed = data.table::fread(highlight.bed)
        colnames(highlight.bed)[1:3] = c("chr", "start", "end")
      }
    }



    ### Build configuration table
    config = data.frame(file = tracks,
                        format = tracks_format,
                        label = track.labels,
                        color = track.colors)




    ########################################################################################################################################################################
    ## Definition of plotting functions for each format of track file
    ########################################################################################################################################################################

    ###### BigWig plot
    plot_bigWig = function(bigWig, genomic.region, color, smooth, span = 0.05, plot.area = T) {

      scores = Rseb::get.single.base.score.bw(region = paste0(genomic.region[1], ":", genomic.region[2], "-", genomic.region[3]),
                                              bigWig = bigWig,
                                              missing.data.as.zero = missing.data.as.zero.bw,
                                              reverse.score = F)
      scores = data.frame(position = 1:length(scores) + as.numeric(genomic.region[2]) - 1, score = scores)  # assign the positions corresponding to the original genomic region


      # plotting of the bigWig
      plot = ggplot(data = scores, aes(x = position, y = score))

      if (smooth == T) {
        if (plot.area == T) {
          plot = plot + geom_ribbon(aes(ymin = 0, ymax = predict(loess(score ~ position, span = span))), fill = color)
        } else {plot = plot + geom_smooth(method = "loess", se = F, formula = y ~ x, span = span, color = color)}
      } else {
        if (plot.area == T) {
          plot = plot + geom_area(show.legend = F, fill = color)
        } else {plot = plot + geom_line(color = color)}
      }

      plot =
        plot +
        scale_x_continuous(expand = c(0,0)) +
        geom_hline(yintercept = 0, color = color) +
        theme_null()

      if (TRUE %in% (scores$score < 0)) { # whe negative values are present, do:
        plot = plot + scale_y_continuous(expand = c(0,0))
      } else {
        plot = plot + scale_y_continuous(expand = c(0,0), limits = c(0, NA))
      }

      return(plot)
    } # end plot_bigWig




      ###### bed plot
      plot_bed = function(bed, genomic.region, color, plot.score = score.bed.shadow) {

        # import bed and coerce to BED6-format
        bed = data.table::fread(bed)

        if (ncol(bed) >= 6) {bed = bed[,1:6]
        } else if (ncol(bed) >= 3 & ncol(bed) < 6) {
          if (ncol(bed) == 3) {
            bed = bed %>% dplyr::mutate(name = paste("region", 1:nrow(bed), sep="_"), score = 0, strand = ".")
          } else if (length(bed) == 4) {
            bed = bed %>% dplyr::mutate(score = 0, strand = ".")
          } else if (length(bed) == 5) {
            bed = bed %>% dplyr::mutate(strand = ".")
          }
        }

        colnames(bed) = c("chr", "start", "end", "name", "score", "strand")


        # select regions overlapping with the genomic.region to plot
        bed = dplyr::filter(.data = bed, chr == genomic.region[1], end > genomic.region[2], start < genomic.region[3]) %>% dplyr::mutate(level = 1)
        bed = Rseb::sort.bed(bed = bed, return.bed = T)

        # Define overlaps
        if (nrow(bed) > 1 & expand.bed == T) {
          level = 1
          continue = TRUE
          bed = dplyr::mutate(.data = bed, overlap = F)

          while (continue) {
            for (i in 2:nrow(bed)) {
              if (bed$level[i] == level) {
                if ((bed$start[i] < bed$end[(i-1)]) & (bed$level[i] == bed$level[i-1])) {
                  bed$overlap[i] = T
                  bed$level[i] = bed$level[i]+1
                }
              }
            }

            continue = (T %in% bed$overlap)
            level = level+1
            bed = dplyr::mutate(.data = bed, overlap = F)
          }
        }

        bed = dplyr::mutate(.data = bed[,1:7], ymin = level-0.2, ymax = level+0.2)


        # generation plot
        if (plot.score == T) {
          plot =
            ggplot(data = bed) +
            geom_rect(data = bed,
                      aes(xmin = start,
                          xmax = end,
                          ymin = ymin,
                          ymax = ymax,
                          alpha = score),
                      fill = color,
                      inherit.aes = FALSE,
                      show.legend = F)
        } else {
          plot =
            ggplot(data = bed) +
            geom_rect(data = bed,
                      aes(xmin = start,
                          xmax = end,
                          ymin = ymin,
                          ymax = ymax),
                      fill = color,
                      inherit.aes = FALSE,
                      show.legend = F)
        }


        plot =
          plot +
          scale_x_continuous(expand = c(0,0), limits = c(as.numeric(genomic.region[2]), as.numeric(genomic.region[3]))) +
          theme_null()

        return(plot)
    }


      # Plot arcs
      plot_arcs =
        function(arcs,
                 arcs.direction,
                 fraction.arc.base,
                 color,
                 genomic.region) {
          # Read arcs
          arcs = data.table::fread(arcs)[,1:6]
          colnames(arcs) = c("chr1", "start1", "end1", "chr2", "start2", "end2")

          arcs =
            arcs %>%
            dplyr::mutate(center1 = (start1 + end1)/2,
                          center2 = (start2 + end2)/2,
                          middle_point = (center1 + center2)/2,
                          r = (((center1 + center2)/2) - center1) * arcs.direction)

          for (i in 1:nrow(arcs)) {if (arcs$center1[i] > arcs$center2[i]) {arcs$r[i] = -arcs$r[i]}}

          plot =
            ggplot(arcs) +
            geom_arc(aes(x0 = middle_point,
                         y0 = 1,
                         r = r,
                         start = pi/2,
                         end = -pi/2),
                     color = color,
                     show.legend = F) +
            geom_rect(data = arcs,
                      aes(xmin = start1,
                          xmax = end1,
                          ymin = max(r) * fraction.arc.base,
                          ymax = max(r) * -fraction.arc.base),
                      fill = color,
                      inherit.aes = FALSE,
                      show.legend = F) +
            geom_rect(data = arcs,
                      aes(xmin = start2,
                          xmax = end2,
                          ymin = max(r) * fraction.arc.base,
                          ymax = max(r) * -fraction.arc.base),
                      fill = color,
                      inherit.aes = FALSE,
                      show.legend = F) +
            scale_x_continuous(expand = c(0,0), limits = c(as.numeric(genomic.region[2]), as.numeric(genomic.region[3]))) +
            scale_y_continuous(expand = c(0,0)) +
            theme_null()

          return(plot)
        }



      # Function to add highlights in the plot
      add.highlight =
        function(p,
                 highlight.bed,
                 highlight.color,
                 highlight.transparency) {
          for (i in c(1:nrow(highlight.bed))) {
            p = p + annotate("rect",
                             xmin = highlight.bed$start[i],
                             xmax = highlight.bed$end[i],
                             ymin = -Inf,
                             ymax = Inf,
                             fill = highlight.color,
                             alpha = highlight.transparency)
          }
          return(p)
        }



      # Plot genomic region annotation (genes)
      gene_annotation_plot =
        autoplot(ensdb,
                 AnnotationFilter::GRangesFilter(GRanges(seqnames = as.numeric(gsub("chr", "", genomic.region[1])),
                                                         IRanges(as.numeric(genomic.region[2]), as.numeric(genomic.region[3])),
                                                         strand = "*")),
                 names.expr = "gene_name",
                 label.color = "#000000", #black
                 color = gene.annotation.color,
                 fill = gene.annotation.color) +
        xlab(paste0(genomic.region[1], ":", genomic.region[2], "-", genomic.region[3])) +
        scale_x_continuous(expand = c(0,0), limits = c(as.numeric(genomic.region[2]), as.numeric(genomic.region[3]))) +
        scale_y_continuous(expand = c(0.1,0)) +
        theme_classic() +
        theme(axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(color = "#000000"))

      if (!is.null(highlight.bed)) {
          gene_annotation_plot = add.highlight(p = gene_annotation_plot, highlight.bed = highlight.bed, highlight.color = highlight.color, highlight.transparency = highlight.transparency)
      }



      ######################################################################################################################################################
      # Generate the single tracks
      track_plot_list = list()

      for (i in 1:nrow(config)) {
        if (config$format[i] == "bigWig") {
          track_plot_list[[i]] = plot_bigWig(bigWig = config$file[i], genomic.region = genomic.region, color = config$color[i], smooth = smooth.bigWig.signal, span = smooth.bigWig.loess.span, plot.area = plot.bigWig.area)
        } else if (config$format[i] == "bed") {
          track_plot_list[[i]] = plot_bed(bed = config$file[i], genomic.region = genomic.region, color = config$color[i])
        } else {
          track_plot_list[[i]] = plot_arcs(arcs = config$file[i], arcs.direction = arcs.direction, fraction.arc.base = fraction.arc.base, color = config$color[i], genomic.region = genomic.region)
        }
      }


      # Add the highlighting if required
      if (!is.null(highlight.bed)) {
        for (i in 1:length(track_plot_list)) {
          track_plot_list[[i]] = add.highlight(p = track_plot_list[[i]], highlight.bed = highlight.bed, highlight.color = highlight.color, highlight.transparency = highlight.transparency)
        }
      }

      track_plot_list[[length(track_plot_list)+1]] = gene_annotation_plot@ggplot


      # Redefine the axis if required by 'grouping'
      if (!is.null(grouping)) {
        for (i in 1:length(grouping)) {
          uniformed_plots = Rseb::uniform.y.axis(plot.list = track_plot_list[grouping[[i]]])
          for (k in 1:length(grouping[[i]])) { track_plot_list[[grouping[[i]][k]]] = uniformed_plots[[k]] }
        }
      }


      # Add axis-range to bigWig tracks
      for (i in 1:length(config$format)) {
        if (config$format[i] == "bigWig") {
          plot_build = ggplot_build(track_plot_list[[i]])

          track_plot_list[[i]] =
            track_plot_list[[i]] +
            annoate("text",
                    label = paste0("[", round(plot_build$layout$panel_params[[1]]$y$limits[1], 2),
                                   " - ", round(plot_build$layout$panel_params[[1]]$y$limits[2], 2), "]"),
                    y = Inf,
                    x = -Inf,
                    color = "#000000",
                    hjust = 0,
                    vjust = 1,
                    size = bigWig.range.label.size)
        }
      }



      # Generate track labels
      track.labels.full = c(config$label,
                            paste0(gsub(".sqlite", "", basename(ensdb@ensdb@dbname)), " (", genome, ")"))
      names(track_plot_list) = track.labels.full

      labels_plot_list = list()

      for (i in 1:length(track_plot_list)) {
        build_plot = ggplot_build(track_plot_list[[i]])
        ymin = build_plot$layout$panel_params[[1]]$y$limits[1]
        ymax = build_plot$layout$panel_params[[1]]$y$limits[2]

        labels_plot_list[[i]] =
          ggplot() +
          annotate("text",
                   label = track.labels.full[i],
                   x = 0,
                   y = mean(c(ymin, ymax)),
                   hjust = 1,
                   vjust = 1,
                   size = track.labels.fontzise,
                   color = "#000000") +
          scale_y_continuous(expand = c(0,0), limits = c(ymin, ymax)) +
          scale_x_continuous(expand = c(0,0), limits = c(-0.1,0)) +
          theme_null()
      }


      #######################################################################################################################
      # Assembling of the multi.plot
      if (is.null(height.ratios)) {
        ratios = c()
        for (i in 1:nrow(config)) {
          if (config$format[i] == "bigWig") {
            ratios[i] = 0.25
          } else if (config$format[i] == "bed") {
            ratios[i] = 0.1
          } else {
            ratios[i] = 0.5}
        }
        ratios = c(ratios, 1) #adding ratio for gene annotation
      } else {ratio = height.ratios}


      multi.track = cowplot::plot_grid(plotlist = Rseb::combine.lists(list(labels_plot_list, track_plot_list)),
                                       ncol = 2,
                                       rel_heights = ratios,
                                       rel_widths = width.ratios,
                                       byrow = F)



      #################################################################################################################
      # Generating output
      return(list(configuration = rbind(config,
                                        data.frame(file = ensdb@ensdb@dbname,
                                                   format = "genome_annotation",
                                                   label = track.labels.full[length(track.labels.full)],
                                                   color = gene.annotation.color)),
                  highlighted.region = highlight.bed,
                  single.track.list = track_plot_list,
                  single.label.plot.list = labels_plot_list,
                  multi.track.plot = multi.track))
  } # END function
