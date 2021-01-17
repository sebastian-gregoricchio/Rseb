![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/Rseb)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/Rseb)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/Rseb)
![license](https://img.shields.io/badge/License-GPLv3-blue.svg)

## Dependencies/Requirements <img src="Rseb_logo.svg" align="right" height = 120/>
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
devtools::install_github("sebastian-gregoricchio/Rseb")
```
<br /> 

## Documentation
With the package a [PDF manual](https://sebastian-gregoricchio.github.io/Rseb/Rseb_0.1.4_manual.pdf) is available.

<br /> 

## Package history and releases
The changeLog could be found [here](https://sebastian-gregoricchio.github.io/Rseb/NEWS).

**Old releases**
* [Rseb v0.1.0](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.0)
* [Rseb v0.1.1](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.1)
* [Rseb v0.1.2](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.2)
* [Rseb v0.1.3](https://github.com/sebastian-gregoricchio/Rseb/releases/tag/0.1.3)


<br /> 
<br /> 

-----------------
## Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

## License
This package is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).

<br /> 

##### Contributors
![[https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true]](https://sebastian-gregoricchio.github.io/)
