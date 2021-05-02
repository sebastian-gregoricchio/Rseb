## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">")
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
#options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  data("RNAseq", package = "Rseb")
#  RNAseq

## ---- echo=FALSE--------------------------------------------------------------
data("RNAseq", package = "Rseb")

set.seed(floor(runif(1, 1, 1000)))
n = floor(runif(n = 10, min = 1, max = nrow(RNAseq)))
knitr::kable(RNAseq[n,], row.names = F, caption = "DESeq2 results table example")

## ----message=FALSE, warning=FALSE---------------------------------------------
require(dplyr)

RNAseq <-
  RNAseq %>%
  mutate(DE.status = Rseb::DE.status(log2FC = RNAseq$log2FC,
                                     p.value.adjusted = RNAseq$padj,
                                     FC_threshold = 2, # Linear value
                                     FC_NoResp_left = 0.9, # Automatically 0.9 <= FC <= 1/0.9)
                                     p.value_threshold = 0.05,
                                     low.FC.status.label = "DOWN",
                                     high.FC.status.label = "UP",
                                     unresponsive.label = "UNRESPONSIVE",
                                     null.label = "NULL"))


## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(RNAseq[n,], row.names = F, caption = "RNA-seq table with differential expression status")

## ----message=FALSE, warning=FALSE---------------------------------------------
RNAseq.summary.table <-
  RNAseq %>%
  group_by(DE.status) %>%
  summarise(N = n()) %>%
  rbind(c("Total", nrow(RNAseq)))


## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(RNAseq.summary.table, row.names = F, caption = "RNA-seq differential expression summary")

## ----message=FALSE, warning=FALSE, fig.cap='MA-plot of RNA-seq data. <br> Correlation between log~2~(mean normalized expression in all samples and replicates) and log~2~(Fold Change expression) among two conditions are reprresented as dots for each single gene. Up-regulated (FC ≥ 2, P~adj~ < 0.05), down-regulated (FC ≤ 0.5, P~adj~ < 0.05), unresponsive (0.9 ≤ FC ≤ 1.1, P~adj~ ≥ 0.05) or null genes are indicated as pink, green, blue or gray points respectively.', fig.align="center"----
require(ggplot2)

MA.plot <-
  ggplot(data = RNAseq,
         aes(x = log2(baseMean),
             y = log2FC,
             col = DE.status)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#F8766D", "gray30", "#00A5CF", "#00BA38")) +
  ggtitle("MA-plot") +
  theme_classic()

MA.plot

## ----message=FALSE, warning=FALSE, fig.cap='Volcano plot of RNA-seq data. <br> Correlation between log~2~(Foldchange expression) among two conditions and -log~10~(p-value~adjusted~) of the comparisons are reprresented as dots for each single gene. Up-regulated (FC ≥ 2, P~adj~ < 0.05), down-regulated (FC ≤ 0.5, P~adj~ < 0.05), unresponsive (0.9 ≤ FC ≤ 1.1, P~adj~ ≥ 0.05) or null genes are indicated as pink, green, blue or gray points respectively.', fig.align="center"----
volcano.plot <-
  Rseb::volcano(log2FC_data = RNAseq$log2FC,
                padj_data = RNAseq$padj,
                FC_t = 2,
                p_t = 0.05,
                FC_unresponsive_left = 0.9,
                left_label = "DOWN",
                unresponsive_label = "UNRESPONSIVE",
                right_label = "UP",
                null_label = "NULL",
                title = "Volcano",
                point_size = 2)

volcano.plot

## ----message=FALSE, warning=FALSE, fig.cap='Volcano plot of RNA-seq data. <br> Correlation between log~2~(Foldchange expression) among two conditions and -log~10~(p-value~adjusted~) of the comparisons are reprresented as dots for each single gene. Up-regulated (FC ≥ 2, P~adj~ < 0.05), down-regulated (FC ≤ 0.5, P~adj~ < 0.05), unresponsive (0.9 ≤ FC ≤ 1.1, P~adj~ ≥ 0.05) or null genes are indicated as pink, green, blue or gray points respectively. Names of Down-regulated genes are labelled.', fig.align="center"----
volcano.plot.with.names <-
  Rseb::volcano(log2FC_data = RNAseq$log2FC,
                padj_data = RNAseq$padj,
                names = RNAseq$geneName,
                FC_t = 2,
                p_t = 0.05,
                FC_unresponsive_left = 0.9,
                left_label = "DOWN",
                unresponsive_label = "UNRESPONSIVE",
                right_label = "UP",
                null_label = "NULL",
                title = "Volcano",
                point_size = 2,
                right_names = T,
                names_size = 3)

volcano.plot.with.names

## ---- echo=FALSE, out.width = '100%', fig.cap = 'Sequencing signal visualization ([Full size image](https://sebastian-gregoricchio.github.io/Rseb/vignettes/images/Rseb_workflow.svg)). <br> **(A)** The function computeMatrix from deeptools used by [`computeMatrix.deeptools`](#matrix), in the "reference-point mode", extends the center of all the regions in one or more datasets (.bed files) and divides each one in *n* bins of a desired size (eg. 20bp). Then the signal of one or more signals tracks (.bw files) is calculated at each bin and the scores values are stored in a matrix. **(B)** The function [`plot.density.profile`](#densityProfile) calculates the mean/median/sum and the variance for each bin of all the regions in a dataset for each track signal provided (left). Then can be generated a plot (right) showing the profile density around the center of the regions. The plot could be generated by sample (a plot per sample with a different curve for each dataset of regions) or by region (a plot per region with a different curve for each sample). **(C)** The function [`plot.density.summary`](#densitySummary) computes the mean/median/sum of the signal over each single region for each sample. Then violins plot can be generated to reppresent the distribution of all the single density values of the regions by dataset or sample (simalary to [`plot.density.profile`](#densityProfile)). The function computes as well the means statistical comparison whose corresponding significance level/p-value can be directly added to the plots. **(D)** The function [`plot.density.differences`](#densityDifferences) computes the difference of mean/median/sum signal for each single region for each sample in each group. Two graphical rappresentation are available. On top, an [area/line plot](#areaPlot) showing the values of the signal difference for each region which are ranked by these values. On the X-axis it is indicated the number of regions with negative and positive difference. On bottom, a [scatter/correlation plot](#correlationPlot) showing the correlation between the mean/median/sum signal of a region in one sample vs an other one. Positive and negative values of the difference in the signal are indicated by different colors. The function infers the significance for the means difference in the groups for all the pair comparisons and returns it as a list of tables for each group.'----
#magick::image_read(rsvg::rsvg_svg("url_image"))
knitr::include_graphics("https://sebastian-gregoricchio.github.io/Rseb/vignettes/images/Rseb_workflow.svg")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  bed <- Rseb::build.bed(chr = paste("chr", c(round(runif(8,1,23)),"X","Y"), sep=""),
#                         start = round(runif(10,1,100000)),
#                         end = round(runif(10,100001,1000000)),
#                         name = paste("peak_", round(runif(10,1,1000)), sep=""),
#                         strand = "+",
#                         bed.file.name = "/path/to/ouput/file.bed",
#                         return.data.frame = T)

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
bed <- Rseb::build.bed(chr = paste("chr", c(round(runif(8,1,23)), "X", "Y"), sep = ""),
                       start = round(runif(10,1,100000)),
                       end = round(runif(10,100001,1000000)),
                       name = paste("peak_", round(runif(10,1,1000)), sep=""),
                       strand = "+",
                       return.data.frame = T)
knitr::kable(bed, row.names = F, caption = "Bed file example generated by `build.bed`")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  computeMatrix.deeptools(
#    mode = "reference-point",
#    scoreFileName = c("path/to/signal_1.bw", "path/to/signal_3.bw"),
#    regionsFileName = c("path/to/regions_1.bed", "path/to/regions_2.bed", "path/to/regions_3.bed"),
#    outFileName = "/path/to/matrix_file.gz",
#    upstream = 2000, #bp
#    downstream = 2000, #bp
#    binSize = 50, #bp
#    referencePoint = "center",
#    smartLabels = T,
#    missingDataAsZero = T,
#    computeMatrix.deeptools.command = "/home/user/anaconda3/bin/computeMatrix",
#    quiet = T)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  matrix <- read.computeMatrix.file(matrix.file = "/path/to/matrix_file.gz")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  data("deeptools.matrix", package = "Rseb")
#  deeptools.matrix$metadata

## ---- echo=FALSE--------------------------------------------------------------
data("deeptools.matrix", package = "Rseb")

set.seed(floor(runif(1, 1, 1000)))
n = floor(runif(n = 10, min = 1, max = nrow(deeptools.matrix$matrix.data)))

knitr::kable(deeptools.matrix$metadata, row.names = F, caption = "Example of `read.computeMatrix.file` result: metadata")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  deeptools.matrix$data.table

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(deeptools.matrix$matrix.data[n,1:10], row.names = F, caption = "Example of `read.computeMatrix.file` result: matrix.data")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  deeptools.matrix$original.file.path
#  > [1] "/user/path/to/matrix.deepTools.gz"

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
data("deeptools.matrix", package = "Rseb")
density.profile.by.group <-
  Rseb::plot.density.profile(
    matrix.file = deeptools.matrix,
    signal.type = "mean",
    error.type = "sem",
    plot.by.group = T,
    missing.data.as.zero = T,
    y.identical.auto = T,
    text.size = 10,
    axis.line.width = 0.25,
    line.width = 0.5,
    plot.error = T,
    write.reference.points = T,
    plot.vertical.lines = F,
    colors = c("steelblue", "mediumseagreen"),
    n.row.multiplot = 2,
    by.row = T,
    y.lab = "Mean ChIP-seq signal \u00b1 SEM")

str(density.profile.by.group, max.level = 1, give.attr = F)

## ----message=FALSE, warning=FALSE, fig.cap='Density profile. <br> Example of density plot profiles of ChIP-seq normalized signal around the center of four different dataset regions for each sample **(A)** or group **(B)**.', fig.width=6, fig.height=10.5, fig.align="center"----
# Load the example matrix list
data("deeptools.matrix", package = "Rseb")

# Generate the density profile by group
density.profile.by.group <-
  Rseb::plot.density.profile(
    matrix.file = deeptools.matrix,
    signal.type = "mean",
    error.type = "sem",
    plot.by.group = T,
    missing.data.as.zero = T,
    y.identical.auto = T,
    text.size = 10,
    axis.line.width = 0.25,
    line.width = 0.5,
    plot.error = T,
    write.reference.points = T,
    plot.vertical.lines = F,
    colors = c("steelblue", "mediumseagreen"),
    n.row.multiplot = 2,
    by.row = T,
    y.lab = "Mean ChIP-seq signal \u00b1 SEM")

# Generate the density profile by sample
density.profile.by.sample <-
  Rseb::plot.density.profile(
    matrix.file = deeptools.matrix,
    signal.type = "mean",
    error.type = "sem",
    plot.by.group = FALSE,
    missing.data.as.zero = TRUE,
    y.identical.auto = TRUE,
    text.size = 10,
    axis.line.width = 0.25,
    line.width = 0.5,
    plot.error = TRUE,
    write.reference.points = TRUE,
    plot.vertical.lines = FALSE,
    colors = rep(c("indianred", "darkgoldenrod2"), 2),
    line.type = rep(c("solid", "dotted"), each = 2),
    n.row.multiplot = 1,
    legend.position = c(0.2, 0.75),
    y.lab = "Mean ChIP-seq signal \u00b1 SEM")


cowplot::plot_grid(density.profile.by.group$multiplot,
                   density.profile.by.sample$multiplot,
                   nrow = 2, labels = "AUTO", rel_heights = c(2, 1))

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
data("deeptools.matrix", package = "Rseb")
summary.plot.by.group <- 
    Rseb::plot.density.summary(
      matrix.file = deeptools.matrix,
      plot.by.group = TRUE,
      missing.data.as.zero = TRUE,
      signal.type = "sum",
      stat.paired = FALSE,
      stat.hide.ns = FALSE, 
      axis.line.width = 0.5,
      mean.symbol.size = 0.2,
      y.identical.auto = TRUE,
      text.size = 10,
      border.width = 0.25,
      subset.range = c(-1000, 1000), #bp from 0 point
      colors = c("steelblue", "mediumseagreen"),
      legend.position = "right",
      n.row.multiplot = 2,
      by.row = TRUE)

str(summary.plot.by.group, max.level = 1, give.attr = F)

## ----message=FALSE, warning=FALSE, fig.cap='Density summary plot. <br> Violin plots indicate the distribution of signal **(A)** for each dataset comparing all the samples in each dataset, or **(B)** for each sample comparing all the datasets. Horizontal bars with stars indicate the means comparison perfomerd by paired Wilcoxon test. P* < 0.05, P** < 0.01, P*** < 0.001, P**** < 0.0001, ns = not significant. ', fig.width=9, fig.height=5, fig.align="center"----
# Load the example matrix list
data("deeptools.matrix", package = "Rseb")

# Generate the density profile by group
summary.plot.by.group <- 
    Rseb::plot.density.summary(
      matrix.file = deeptools.matrix,
      plot.by.group = TRUE,
      missing.data.as.zero = TRUE,
      signal.type = "sum",
      stat.paired = FALSE,
      stat.hide.ns = FALSE, 
      axis.line.width = 0.5,
      mean.symbol.size = 0.2,
      y.identical.auto = TRUE,
      text.size = 10,
      border.width = 0.25,
      subset.range = c(-1000, 1000), #bp from 0 point
      colors = c("steelblue", "mediumseagreen"),
      legend.position = "right",
      n.row.multiplot = 2,
      by.row = TRUE)

# Generate the density profile by sample
summary.plot.by.sample <- 
    Rseb::plot.density.summary(
      matrix.file = deeptools.matrix,
      plot.by.group = FALSE,
      missing.data.as.zero = TRUE,
      signal.type = "sum",
      stat.paired = FALSE,
      stat.hide.ns = FALSE, 
      axis.line.width = 0.5,
      mean.symbol.size = 0.2,
      y.identical.auto = TRUE,
      text.size = 10,
      border.width = 0.25,
      subset.range = c(-1000, 1000), #bp from 0 point
      colors = c("steelblue", "mediumseagreen"),
      legend.position = "right",
      n.row.multiplot = 2,
      by.row = TRUE)

cowplot::plot_grid(summary.plot.by.group$multiplot,
                   summary.plot.by.sample$multiplot,
                   nrow = 1, labels = "AUTO", rel_widths = c(3, 2))

## ----message=FALSE, warning=FALSE, fig.cap='Combined density summary plot. <br> Violin plots indicate the distribution of signal for all datasetes for all the samples comparing **(A)** the datasetes for each sample or **(B)** the samples for each dataset.', fig.width=8.5, fig.height=6, fig.align="center"----
cowplot::plot_grid(summary.plot.by.group$summary.plot.samples,
                   summary.plot.by.group$summary.plot.regions,
                   nrow = 2, labels = "AUTO")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  summary.plot.by.sample$means.comparisons

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(summary.plot.by.sample$means.comparisons, row.names = F, caption = "Example of comparisons table for the plots by sample")

## ----message=FALSE, warning=FALSE, fig.cap='Correlation plot. <br> Correlation of the signal over each region in each datasetes for one sample is correlated with the signal in the same region in the second sample. Blue dots correspond to regions with a signal higher in Sample1 while green dots indicate regions with a higher signal in Sample2, while regions whose signal is unchanged among the two samples are indicated in violet. The diagonal dashed line corresponds to the $y = x$ equation. Scores of each dot are projectied perpendicularly to each axis and indicated by a green/blue line. Orange line indicates the correlation function that exists between the two samples with the relative SEM. Correlation equation and R-squared for the regression are indicated on the plot.', fig.width=7, fig.height=7, fig.align="center"----
# Load the example matrix list
data("deeptools.matrix", package = "Rseb")

# Generate the correlation.plot
correlation.plot <- 
  Rseb::plot.density.differences(matrix.file = deeptools.matrix,
                                 inverted.comparisons = T,
                                 missing.data.as.zero = T,
                                 signal.type = "sum",
                                 stat.paired = T,
                                 points.size = 0.25,
                                 axis.line.width = 0.25,
                                 area.y.identical.auto = T,
                                 text.size = 10,
                                 correlation.show.equation = T,
                                 correlation.correlation.line.width = 0.5,
                                 correlation.correlation.line.color = "#ED8141",
                                 subset.range = c(-1000, 1000), #bp from 0 point
                                 colors = c("steelblue", "mediumseagreen", "purple"),
                                 legend.position = c(0.8, 0.25),
                                 error.type = "sem")

cowplot::plot_grid(plotlist = Rseb::uniform.y.axis(Rseb::uniform.x.axis(correlation.plot$correlation.plot.byGroup.list)),
                   nrow = 2, byrow = T)

## ----message=FALSE, warning=FALSE, fig.cap='Area plot. <br> Distribution of the difference between Sample1 (blue) and Sample2 (green) of the signals within each region for each dataset. The dots represent the signal difference value over each region. Blue and green dots signal intensity higher in Sample1 and Sample2 respectively. Regions with unchanged signal among the two samples are indicated in purple.', fig.width=7, fig.height=6, fig.align="center"----
# Load the example matrix list
data("deeptools.matrix", package = "Rseb")

# Generate the area.plot
areaDifference.plot <- 
  Rseb::plot.density.differences(matrix.file = deeptools.matrix,
                                 inverted.comparisons = T,
                                 missing.data.as.zero = T,
                                 signal.type = "sum",
                                 stat.paired = T,
                                 points.size = 0.25,
                                 axis.line.width = 0.25,
                                 area.y.identical.auto = T,
                                 text.size = 10,
                                 subset.range = c(-1000, 1000), #bp from 0 point
                                 colors = c("steelblue", "mediumseagreen", "purple"),
                                 legend.position = c(0.2, 0.80),
                                 error.type = "sem")

cowplot::plot_grid(plotlist = Rseb::uniform.y.axis(areaDifference.plot$area.plot.byGroup.list),
                   nrow = 2, byrow = T)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  # Define the gene list
#  gene_list <- c("Spi1", "Idh1", "Bcl2l11", "Mcl1", "Polr2a", "Hdac1")
#  
#  # Generate the script
#  Rseb::IGVsnap(loci_vector = gene_list,
#                input_type = "genes",
#                biomart = "ensembl",
#                dataset = "mmusculus_gene_ensembl",
#                reference_genome = "mm10",
#                fivePrime = 1000,
#                threePrime = 1000,
#                snap_names = gene_list,
#                IGV_batch_file = "path/to/output/script.txt",
#                snap_directory = "path/to/directory/in/wich/the/images/will/be/generated",
#                snap_image_format = "png",
#                maxPanelHeight = 1500)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  Rseb::uniform.x.axis(plot.list = list.of.plots,
#                       x.min = 0,
#                       x.max = NA,
#                       ticks.each = 0.5,
#                       digits = 1)

## ---- eval=FALSE--------------------------------------------------------------
#  # Loading the table
#  data("CNV.data", package = "Rseb")
#  CNV.data

## ---- echo=FALSE--------------------------------------------------------------
data("CNV.data", package = "Rseb")

knitr::kable(CNV.data[round(runif(n = 10, min = 1, max = nrow(CNV.data))),], row.names = F, caption = "Example of CNV annotation table")

## ----eval=FALSE---------------------------------------------------------------
#  CNV.gain.loss <- dplyr::filter(.data = CNV.data,
#                                 Rseb::grepl.data.frame(data.frame = CNV.data,
#                                                        pattern = "gain/loss"))
#  print(CNV.gain.loss)

## ---- echo=FALSE--------------------------------------------------------------
CNV.gain.loss = dplyr::filter(.data = CNV.data, Rseb::grepl.data.frame(data.frame = CNV.data, pattern = "gain/loss"))

knitr::kable(CNV.gain.loss, row.names = F, caption = 'CNV annotation table filtered for the pattern "gain/loss"')

