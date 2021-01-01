# opcr
Process and Plot Data from an Optical Plankton Counter (OPC)

## Installation

You can install the latest version of opcr from github:
``` r
library(devtools)
install_github('hansenjohnson/opcr')
library(opcr)
```

## Example

Here's a quick example of how to process and plot OPC data. See the `intro` vignette for more details.

``` r
library(opcr)

# process opc data file (requires raw data file)
opc = opc_process('data/raw/OPC001.D00')

# plot abundance histogram
opc_plot_histogram(opc)
```
