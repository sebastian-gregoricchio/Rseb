# Rseb <img src="https://sebastian-gregoricchio.github.io/Rseb/Rseb_logo.svg" align="right" height = 150/>

An R-package for daily tasks required to handle biological data as well as avoid re-coding of small functions for quick but necessary data managing.

## Dependencies/Requirements
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

### deepTools
Certain functions of this package require that `deeptools` is installed on your system. For more information see the [deepTools](https://deeptools.readthedocs.io/en/develop/content/installation.html) installation page.
* **Installation via `conda`**: `conda install -c bioconda deeptools`
* **Command line installation using `pip`**: `pip install deeptools`, All python requirements should be automatically installed.


### bedTools
Certain functions of this package require that `bedtools` is installed on your system. For more information see the [bedTools](https://bedtools.readthedocs.io/en/latest/content/installation.html) installation page.
* **Installation via `conda`**: `conda install -c bioconda bedtools`
* **Command line installation**:
    - Fedora/Centos: `yum install BEDTools`
    - Debian/Ubuntu: `apt-get install bedtools`
    - Homebrew (MacOS): `brew tap homebrew/science; brew install bedtools`
    - MacPorts: `port install bedtools`


## Installation
```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
## install.packages("devtools")
## devtools::install_github("r-lib/devtools")

# Install the Rseb package
devtools::install_github("sebastian-gregoricchio/Rseb")
```

## Documentation
With the package a [PDF manual](https://sebastian-gregoricchio.github.io/Rseb/Rseb_0.1.1_manual.pdf) is available.


## Package history and releases
The changeLog could be found [here](https://github.com/sebastian-gregoricchio/Rseb/blob/main/NEWS.md).

**Old releases**
* [Rseb v0.1.0](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.0)



-----------------
## Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

## License
This package is under a [GNU General Public License (version 3)](https://github.com/sebastian-gregoricchio/Rseb/blob/main/LICENSE.md/LICENSE.md).
