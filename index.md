![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/Rseb)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/Rseb)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/Rseb)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/Rseb?style=social)](https://github.com/sebastian-gregoricchio/Rseb/fork)
<!---![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/Rseb/total.svg)--->

## Introduction [<img src="https://sebastian-gregoricchio.github.io/Rseb/Rseb_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/Rseb)

The concept behind the `Rseb` (*R*-package for *S*implified *E*nd-to-end data *B*uild-up) is to provide a toolkit that allows the automation of different type of tasks avoiding retyping of code and loss of time. Furthermore, the advantage is that, in most of the cases, the functions are built in R-language making it suitable for all the operating systems

From a more biological point of view, this package simplifies many downstream analyses of high-throughput data that otherwise would require many hours of code typing and case-to-case adaptation. Moreover, most of the functions aimed to visualize these kind of data are thought to provide a high level of possible customization with a large number of graphical parameters compared to the commonly used already available tools. Another advantage of this package is that it offers multiple methods, with a corresponding visualization, to quantify the difference of signal between samples, in a qualitatively and/or quantitatively way, without any additional coding.

The guide is divided in three parts: the first one will explore some analyses and visualization of RNA-seq data, the second one the representation and quantification of targeted sequencing data (ChIP-seq, ATAC-seq, etc.), while the last part is focused on some of the "general" tools available.


### Citation
If you use this package, please cite:

<div class="warning" style='padding:2.5%; background-color:#ffffee; color:#787878; margin-left:5%; margin-right:5%; border-radius:15px;'>
<span>
<font size="-0.5">

<div style="margin-left:2%; margin-right:2%; text-align: justify">
... the publication is under revision ...
</div>
</font>

</span>
</div>

<br>


## Dependencies/Requirements
### Bioconductor libraries
Some functions of this package require `Bioconductor` libraries. These functions should install automatically with the package.
However, if you prefere to manually install `Bioconductor` and required packages proceede with the `Bioconductor` installation:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

The required Biocondutor packages are: `Biostrings`, `biomaRt`, `diffloop`, `GO.db`, `rtracklayer`, `GenomicRanges`, `AnnotationFilter`, `EnsDb.Hsapiens.v75`, `EnsDb.Hsapiens.v86`, `EnsDb.Mmusculus.v79`.
To install it manually proceed with:

```r
BiocManager::install(c("Biostrings", "biomaRt", "diffloop", "GO.db", "rtracklayer", "GenomicRanges", "AnnotationFilter", "EnsDb.Hsapiens.v75", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79"))
```
<br />

#### deepTools
Certain functions of this package require that `deeptools` is installed on your system. For more information see the [deepTools](https://deeptools.readthedocs.io/en/develop/content/installation.html) installation page.
* **Installation via `conda`**: `conda install -c bioconda deeptools`
* **Command line installation using `pip`**: `pip install deeptools`, All python requirements should be automatically installed.

<br />

#### bedTools
Certain functions of this package require that `bedtools` is installed on your system. For more information see the [bedTools](https://bedtools.readthedocs.io/en/latest/content/installation.html) installation page.
* **Installation via `conda`**: `conda install -c bioconda bedtools`
* **Command line installation**:
    - Fedora/Centos: `yum install BEDTools`
    - Debian/Ubuntu: `apt-get install bedtools`
    - Homebrew (MacOS): `brew tap homebrew/science; brew install bedtools`
    - MacPorts: `port install bedtools`

<br />

## Installation
```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
## install.packages("devtools")
## devtools::install_github("r-lib/devtools")

# Install the Rseb package
devtools::install_github("sebastian-gregoricchio/Rseb",
			 build_manual = TRUE,
                         build_vignettes = TRUE)
```
<br />

## Documentation
With the package a [PDF manual](https://sebastian-gregoricchio.github.io/Rseb/Rseb_manual.pdf) and a [vignette](https://sebastian-gregoricchio.github.io/Rseb/doc/Rseb.overview.vignette.html) are available.


<br />

## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/Rseb/NEWS).

<br />
<br />

-----------------
## Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

## License
This package is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).


<br />

##### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/Rseb/NEWS).

<br />
<br />

-----------------
## Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

## License
This package is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).


<br />

##### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
