#' @title GSEA plotter
#'
#' @description Function to plot GSEA results (see \href{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}{clusterprofiler}).
#'
#' @param gsea.results A \code{gseaResult} object as generate by \href{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}{clusterprofiler}.
#' @param geneset.id Numeric value or a string identifying the Nth geneSet (numeric) or a specific id (string) if geneSet in the result table. Default value: \code{NULL}, which returns the ordered list of available geneSets.
#' @param enrichment.geom String indicating the type of graph to use to plot the enrichment scores. Possible options: 'line', 'lines', 'dot', 'dots', 'point', 'points' (case insensitive). Default: \code{'line'}.
#' @param geneset.colors String vector indicating a list of any R-supported color to be used for the enrichment score plot per each dataset. Default: \code{'NULL'} (automatic rainbow colors).
#' @param enrichment.geom.size Numeric value indicating the size of the line, or dots, used in the enrichment score plot. Default: \code{1}.
#' @param enrichment.plot.zero.line Logical value to indicated whether to plot an horizontal line at 0 in the enrichment score plot. Default: \code{FALSE}.
#' @param enrichment.zero.line.color String indicating any R-supported color to be used for the 0-line in the enrichment score plot (active when \code{enrichment.plot.zero.line = TRUE}). Default: \code{'gray'}.
#' @param enrichment.zero.line.width Numeric value indicating the line width of the 0-line in the enrichment score plot (active when \code{enrichment.plot.zero.line = TRUE}). Default: \code{0.5}.
#' @param geneset.segments.width Numeric value indicating the line width of the geneSet vertical segments. Default: \code{0.3}.
#' @param geneset.segments.color String indicating any R-supported color to be used for the geneSet segments. Default: \code{'black'}.
#' @param ranking.color String indicating any R-supported color to be used for the ranked list plot (histogram). Default: \code{'gray'}.
#' @param gradient.colors Two-values string vector indicating the shadows of palettes to use for the genset gradient. Possible values: 'Blues', 'Greens', 'Greys', 'Oranges', 'Purples', 'Reds'. Default: \code{c('Reds', 'Blues')}.
#' @param table.font.size Numeric value to indicate the font size to use in the stat table plot. Default: \code{8}.
#' @param title.position String indicating the position of the title: 'left', 'center', 'right'. Default: \code{'center'}.
#' @param title String indicating the title to use. Default: \code{NULL} (no title).
#' @param combined.plot.height.rations Numeric vector with 4 values used to define the 'real heights' of each panel. The order corresponds to: enrichment score panel, geneset panel (all together), rank panel, stat table panel. Default \code{c(1, 0.4, 0.6, 0.5)}.
#' @param image.file.name String indicating the full path for the export of a pdf file of the combined plot. Default: \code{NULL}, no plot will be exported.
#' @param image.width Numeric value to indicate the width (in inches) to use for the exported pdf file. Active only when \code{image.file.name} is not \code{NULL}. Default: \code{9}.
#' @param image.height Numeric value to indicate the height (in inches) t use for the exported pdf file. Active only when \code{image.file.name} is not \code{NULL}. Default: \code{6}.
#' @param return.all.objects Logical value to indicate whether the function should return only the combined plot (ggplot object), or all the different panels and the combined plot in a list. Default: \code{FALSE} (only combined plot).
#'
#' @return Either a ggplot-object with the final combined plot, or a list with the three panels separated and the combined plot: \code{list(enrichment.panel, geneset.panel.list - list with one element per geneset -, rank.panel, stat.table, combined.plot)}.
#'
#'
#' @examples
#' data(geneList, package = "DOSE")
#'
#' msigdb_hallmarks =
#'   msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
#'   dplyr::select(gs_name, human_entrez_gene)
#'
#' gsea_H = clusterProfiler::GSEA(geneList = geneList,
#'                                TERM2GENE = msigdb_hallmarks,
#'                                minGSSize = 3,
#'                                maxGSSize = 800,
#'                                pAdjustMethod = "BH",
#'                                pvalueCutoff = 0.05,
#'                                verbose = T)
#'
#' plot.multi.gsea(gsea_H, geneset.id = c("HALLMARK_ADIPOGENESIS", "HALLMARK_CELL_CYCLE"))
#' plot.multi.gsea(gsea_H, geneset.id = c(1,2,7))
#'
#' @export plot.multi.gsea




plot.multi.gsea =
  function(gsea.results,
           geneset.id.list = NULL,
           enrichment.geom = "line",
           geneset.colors = NULL,
           enrichment.geom.size = 1,
           enrichment.plot.zero.line = FALSE,
           enrichment.zero.line.color = "gray",
           enrichment.zero.line.width = 0.5,
           geneset.segments.width = 0.3,
           geneset.segments.color = "black",
           ranking.color = "gray",
           gradient.colors = c("Reds", "Blues"),
           table.font.size = 8,
           title.position = "center",
           title = NULL,
           combined.plot.height.rations = c(1, 0.4, 0.6, 0.5),
           image.file.name = NULL,
           image.width = 9,
           image.height = 6,
           return.all.objects = FALSE) {


    #----------------------------------------#
    # Check if Rseb is up-to-date            #
    Rseb::actualize(update = F, verbose = F) #
    #----------------------------------------#


    # Libraries
    require(ggplot2)


    # collect results
    results = data.frame(gsea.results@result)


    # check geneset.id.list
    check.id =
      function(geneset.id){
        if (is.numeric(geneset.id)) {
          if (geneset.id <= nrow(results) & geneset.id > 0) {gene.set = results$ID[geneset.id]
          } else {return("absent")}
        } else if (is.character(geneset.id)) {
          if (geneset.id %in% results$ID) {gene.set = geneset.id
          } else {return("absent")}
        } else {return("absent")}
        return(gene.set)
      }

    gene.set.list = sapply(geneset.id.list, function(x){check.id(x)}, USE.NAMES = F)

    if ("absent" %in% gene.set.list) {
      warning("The genset.id must be a string or a numeric value among the ones above.")
      print(data.frame(gene.set.ID = results$ID))
      return(invisible())
    }


    # Check graphic parameters
    if (!(tolower(enrichment.geom) %in% c("line", "lines", "dot", "dots", "point", "points"))) {
      return(warning("The 'enrichment.geom' must be one among: 'line', 'lines', 'dot', 'dots', 'point', 'points'."))
    }


    if (tolower(title.position) == "right") {
      title.position = 1
    } else if (tolower(title.position) == "left") {
      title.position = 0
    } else {
      title.position = 0.5
    }


    if (!is.null(geneset.colors)) {
      if (length(geneset.colors) >= length(gene.set.list)) {
        if (!is.null(names(geneset.colors))) {
          set.colors = geneset.colors[1:length(gene.set.list)]
        } else {
          set.colors = geneset.colors[geneset.id.list]
        }
      } else {
        set.colors = rainbow(length(gene.set.list))
      }
    } else {
      set.colors = rainbow(length(gene.set.list))
    }

    names(set.colors) = gene.set.list



    ### Generate a list of plots and get the single plots
    all.enrichment.panels = list()
    all.set.panels = list()

    for (i in 1:length(gene.set.list)) {
      gsea = Rseb::plot.gsea(gsea.results = gsea.results,
                             geneset.id = gene.set.list[i],
                             enrichment.geom = enrichment.geom,
                             enrichment.color = set.colors[i],
                             geneset.segments.color = set.colors[i],
                             gradient.colors = gradient.colors,
                             geneset.segments.width = geneset.segments.width,
                             title = NA,
                             return.all.objects = T)


      all.enrichment.panels[[i]] = ggplot_build(gsea$enrichment.panel)$data[[1]]
      all.enrichment.panels[[i]]$`Geneset ID` = gene.set.list[i]

      all.set.panels[[i]] = gsea$geneset.panel

      if (i == 1) {rank.panel = gsea$rank.panel}
    }

    names(all.set.panels) = gene.set.list

    enrichment.panel.tb = do.call(rbind, all.enrichment.panels)


    # Replot enrichments
    if (enrichment.geom %in% c("line", "lines")) {
      enrichment.panel =
        ggplot(data = enrichment.panel.tb,
               aes(x = x,
                   y = y,
                   color = `Geneset ID`)) +
        geom_line(linewidth = enrichment.geom.size) +
        scale_color_manual(values = set.colors)
    } else {
      enrichment.panel =
        ggplot(data = enrichment.panel.tb,
               aes(x = x,
                   y = y,
                   color = `Geneset ID`)) +
        geom_point(size = enrichment.geom.size) +
        scale_color_manual(values = set.colors)
    }

    enrichment.panel =
      enrichment.panel +
      ggtitle(title) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      geom_vline(xintercept = results$rank[which(results$ID %in% gene.set.list)],
                 linetype = "dotted",
                 color = set.colors) +
      xlab(NULL) +
      ylab("Running Enrichment Score (ES)") +
      theme_classic() +
      theme(axis.text.y = element_text(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(c(0,0,0,0)),
            plot.title = element_text(hjust = title.position))

    if (enrichment.plot.zero.line == TRUE) {
      enrichment.panel =
        enrichment.panel +
        geom_hline(yintercept = 0,
                   linetype = 1,
                   linewidth = enrichment.zero.line.width,
                   color = enrichment.zero.line.color)
    }


    ### Add max.rank lines to the rank panel
    rank.panel =
      rank.panel +
      geom_vline(xintercept = results$rank[which(results$ID %in% gene.set.list[-1])],
                 linetype = "dotted",
                 color = set.colors[-1])

    rank.panel$layers[[4]] = NULL #removing rank annotation



    ### Create a stat table
    stat = results[which(results$ID %in% gene.set.list),c(1,3,4,5,6,7,8,9)]
    stat$pvalue = signif(stat$pvalue, 3)
    stat$p.adjust = signif(stat$p.adjust, 3)
    stat[,7] = signif(stat[,7], 3)
    rownames(stat) = which(results$ID %in% gene.set.list)

    stat.tb.plot =
      ggplot() +
      annotation_custom(gridExtra::tableGrob(stat, theme = gridExtra::ttheme_minimal(base_size = table.font.size)), ) +
      theme_void()



    # COMBINING the plots
    combined.plot =
      cowplot::plot_grid(plotlist = Rseb::combine.lists(list(list(enrichment.panel), all.set.panels, list(rank.panel, stat.tb.plot))),
                         ncol = 1,
                         align = "v",
                         axis = "tblr",
                         rel_heights = c(combined.plot.height.rations[1],
                                         rep(combined.plot.height.rations[2]/(length(all.set.panels)), length(all.set.panels)),
                                         combined.plot.height.rations[3],
                                         combined.plot.height.rations[4]))


    # Export image if required
    if (!is.null(image.file.name)) {
      pdf(file = image.file.name, width = image.width, height = image.height)
      print(combined.plot)
      invisible(dev.off())
    }


    if (return.all.objects == T) {
      return(list(enrichment.panel = enrichment.panel,
                  geneset.panel.list = all.set.panels,
                  rank.panel = rank.panel,
                  stat.table = stat.tb.plot,
                  combined.plot = combined.plot))
    } else {
      return(combined.plot)
    }

  } # END FUNCTION
