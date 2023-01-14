---
title: "changeLog"
---

# Rseb [<img src="https://sebastian-gregoricchio.github.io/Rseb/Rseb_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/Rseb)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/Rseb)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/Rseb?style=social)](https://github.com/sebastian-gregoricchio/Rseb/fork)


#### [v0.3.2](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.3.2) - January 14<sup>th</sup> 2023
* Update of `evaluate.heterogeneity` to enable the log1p signal in the heatmap and avoid scientific annotation of the rank numbers.

<br />

#### [v0.3.1](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.3.1) - December 10<sup>th</sup> 2022
* New functions added: `evaluate.heterogeneity`, `incucyte`, `read.computeMatrix.file`, `intersect.regions`, `venn.overlap` and `closest.regions`
* Updated the dependencies to be compatible with new functions: *scales*, *grDevices*, *R.utils*, *eulerr*, *IRanges* and *pryr*
* Added the possibility to export a file with the session info in the function `pkg.version`
* Bug fixing for *abs(number) < 1* in `floating.floor` and `floating.ceiling`
* Optimized alignment of the plot for `plot.density.differences`, `plot.density.profile`, `plot.density.profile.smooth` and `plot.density.summary`

<br />

#### [v0.3.0](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.3.0) - August 15<sup>th</sup> 2022
* Improvement of `sort.bed` usage: a) removed redundant "export.bed" option (now it is sufficient to add a file name in the "export.file.name" parameter if exporting is wished), b) added the "unique.regions" option, c) minor bug fixed in `sort.bed` for input class check
* Added option 'verbose' for `sort.bed` and `collapse.bed` functions
* Added citation file to cite the article associated to this package
* Improvement of the vignette
* Changed the default value for command string parameter for `intersect.bed` and `computeMatrix.deeptools` functions
* Bug-fixing for the handling of data.frames inputs in the function `density.matrix`
* Bug-fixing for the handling of bigWig inputs that do not contain 'chr' in the seqnames/chromosome names for the function `density.matrix`
* Added *stranded* option in the function `density.matrix`, allowing the possibility to reverse the signal of strand '-' regions
* Addition of *diffloop* dependence used to handle bigWig inputs that do not contain 'chr' in the seqnames for the function `density.matrix`
* Added new function `get.single.base.score.bw`
* Added new function `plot.density.profile.smooth`
* Added new function `genomic.tracks`
* Added new function `qPCR.rna.exp`
* Added new function `qPCR.results` example
* Updated dependencies to satisfy `genomic.tracks` and `qPCR.rna.exp` requirements: *colorspace*, *ggbio*, *ggforce*, *GenomicRanges*, *AnnotationFilter*, *EnsDb.Hsapiens.v75*, *EnsDb.Hsapiens.v86*, *EnsDb.Mmusculus.v79*, *openssl*, *readxl*, *jpeg*

<br />

#### [v0.2.0](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.2.0) - April 13<sup>th</sup> 2021
* Updated information in DESCRIPTION file (new suggested: `ggpmisc`, `prettydoc`, `knitr`, `rmarkdown`, `stats`)
* Added [vignette](https://sebastian-gregoricchio.github.io/Rseb/doc/Rseb.overview.vignette.html) to describe the main functions of the package
* Added function `density.matrix` (alternative to `computeMatrix.deeptools`)
* Added example datasets `RNAseq`, `deeptools.matrix`, `CNV.table`
* New function `plot.density.differences`
* Updtate `actualize` function to build vignettes and bugs fixed
* Optimized the `plot.density.summary` and `plot.density.profile` functions in order to not interfer with other packages
* Optimized the `plot.density.summary` and `plot.density.profile` functions in order to keep the input order of the samples and groups in the plots colors/legend
* Updated `plot.density.profile` to change legend title when 'plot.by.group==FALSE'
* Updated `plot.density.summary` to compute the means comparisons even when not required to show on the plot and indicate whether the statistical test performed was paired or not
* Bug fixed in matrix loading for `plot.density.summary` function
* Better handling of Rseb version visualization in `pkg.version` output

<br />

#### [v0.1.7](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.7) - February 20<sup>th</sup> 2021
* Updated information in DESCRIPTION file
* Updated the dependencies to require an R version >= 4.0.0
* Moved all "depending" packages in "imports" to avoid attaching of all packages. Morover, added `ggpubr` and removed `labeling`
* Bug fixing in `actualize` and `uniform.x.axis` functions
* Optimization for better definition of the breaks in the functions `uniform.x.axis` and `uniform.y.axis`
* Optimization of package handling for `plot.density.summary` and `plot.density.profile`
* Optimized color handling for summary plots in `plot.density.summary` function
* Now computation and plot of the statistical mean comparisons is available for `plot.density.summary`
* In the `plot.density.summary` function a density profile corresponding to the summary plot subset range is generated
* The data.table of `plot.density.summary` now contains also the orignal positions used to generate the matrix
* Optimization of `IGVsnap` function and improvment of relative manual information
* `IGVsnap` allows the addiction of a delay time between snap generation

<br />

#### [v0.1.6](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.6) - February 10<sup>th</sup> 2021
* Added the function `actualize` and integrated in each function to check everytime if `Rseb` package is up-to-date
* New functions `plot.density.summary`, `floating.ceiling`, `floating.floor`
* Updated the dependencies to add `rvcheck` and `curl`
* Optimization for different labeling and better handle of Y-axis autoscale in `plot.density.profile` function
* Optimization of `uniform.x.axis` and `uniform.y.axis` with the new functions `floating.ceiling` and `floating.floor` and to support a given number of digits

<br />

#### [v0.1.5](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.5) - January 29<sup>th</sup> 2021
* `DE.status` now supports also NAs for the FoldChange
* Added functions `color.gradient` and `is.color`
* Bugs fixed in `volcano`

<br />

#### [v0.1.4](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.4) - January 15<sup>th</sup> 2021
* Bugs fixed in `plot.density.profile`
* Bugs fixed in `build.bed`
* Bugs fixed in `intersect.bedtools`
* Updated dependencies required (`labeling`)
* Added functions `uniform.x.axis` and `uniform.y.axis`

<br />

#### [v0.1.3](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.3) - January 11<sup>th</sup> 2021
* Optimization of `volcano` function: parameter 'font_size' added
* Bugs fixed in `computeMatrix.deeptools`
* Added function `build.bed`
* Added R version dependence (set to R >= 3.2.0) to be compatible with `dplyr`

<br />

#### [v0.1.2](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.2) - January 7<sup>th</sup> 2021
* Optimization of `IGVsnap` function
* Added the function `pkg.check`

<br />

#### [v0.1.1](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.1) - January 6<sup>th</sup> 2021
* Some bug-fixing
* Added the function `read.computeMatrix.file`
* Optimized function `plot.density.profile` to be compatibile with outputs from `read.computeMatrix.file`

<br />

#### [v0.1.0](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.0) - January 5<sup>th</sup> 2021
First release.




<br />
<br />

-----------------------------------------------------------------------

##### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
