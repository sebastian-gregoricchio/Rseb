# Dependencies/Requirements <img src="Rseb_logo.svg" align="right" height = 120/>
Some functions of this package require `Bioconductor` libraries. These functions should install automatically with the package. Anyway, install `Bioconductor` repository is recommended. 

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

The possibile required packages are: `Biostrings`, `biomaRt`, `GO.db`, `rtracklayer`.
To install them directly you need to have previously installed `BiocManager` and then:

```r
BiocManager::install(c("Biostrings", "biomaRt", "GO.db", "rtracklayer"))
```

# Installation
```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
## install.packages("devtools")
## devtools::install_github("r-lib/devtools")

# Install the Rseb package
devtools::install_github("sebastian-gregoricchio/Rseb")
```

# Documentation
With the package a [PDF manual](https://github.com/sebastian-gregoricchio/Rseb/blob/main/Rseb_0.1.0_manual.pdf) is available.

# Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

# License
This package is under a [GNU General Public License (version 3)](https://github.com/sebastian-gregoricchio/Rseb/blob/main/LICENSE.md/LICENSE.md).