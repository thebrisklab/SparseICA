# SparseICA

Sparse ICA (Sparse Independent Component Analysis) is a novel ICA method that enables sparse estimation of independent source components.

## Installation

We assume you are running R 4.1.0 or newer. There is no guarantee for backward or forward comparability. Please raise the issue on GitHub if something breaks.

The following R packages are required:

- Rcpp (>= 1.0.9)
- RcppArmadillo (>= 0.12.0.0.0)
- MASS (>= 7.3-58.1)
- irlba (>= 2.3.5)
- clue (>= 0.3)
- devtools

You can install them by running this code:

```r
if(!require(c("Rcpp","RcppArmadillo","MASS","irlba","clue","devtools"))){
    install.packages(c("Rcpp","RcppArmadillo","MASS","irlba","clue","devtools"))
}
```

Then you can install Sparse ICA from github with:

``` r
library(devtools)
install_github("thebrisklab/SparseICA")
```

## Quick start guide

### Load package and example data

``` r
# Load the package
library(SparseICA)

# Load the example data
data(example_sim123)
```

### Visualization of example data

The data has




