# opcr
Process and Plot Data from an Optical Plankton Counter (OPC)

## Overview

An Optical Plankton Counter (OPC) is a tried and true piece of oceanographic equipment designed to count and size small particles in the water. The raw data are stored in a binary file format (typically with a .D00 file extension) designed by the manufacturer, Focal Technologies. This package provides utilities for reading in D00 files, converting them to a nice R format, and deriving and plotting relevant metrics.

## Installation and use

You can install the latest version of opcr from github as follows. I recommend building the vignettes as well so you can see some quick examples of how to use the package.

``` r
# install the latest version from github
devtools::install_github('hansenjohnson/opcr', build_vignettes = TRUE)
library(opcr)

# check out the vignette for simple examples
vignette('intro',package='opcr')
```
