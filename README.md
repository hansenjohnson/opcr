# opcr
Process and plot data from an Optical Plankton Counter (OPC)

[![DOI](https://zenodo.org/badge/326001520.svg)](https://zenodo.org/badge/latestdoi/326001520)

## Overview

An Optical Plankton Counter (OPC) is an oceanographic instrument designed to count and size small particles in the water. The raw data are stored in a binary file format (typically with a .D00 file extension) designed by the manufacturer, Focal Technologies. This package provides utilities for reading in D00 files, converting them to a nice R format, and computing and plotting relevant metrics.

## Installation and use

You can install the latest version of opcr from github as follows. I recommend building the vignettes as well so you can see some quick examples of how to use the package.

``` r
# install the latest version from github
devtools::install_github('hansenjohnson/opcr', build_vignettes = TRUE)
library(opcr)

# check out the vignette for simple examples
vignette('opcr')
```

## Citation

Please cite as follows (update year and version as needed):

> Johnson, HD (2021). opcr: process and plot data from an optical plankton counter, v1.0.0. Zenodo. DOI: 10.5281/zenodo.5760977
