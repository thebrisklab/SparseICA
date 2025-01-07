SparseICA: Sparse Independent Component Analysis <img src="fig/sticker.png" width="120" align="right" />
===================================================

<br>

Sparse ICA (Sparse Independent Component Analysis) is a novel ICA method that enables sparse estimation of independent source components.

## Installation

We assume you are running R 4.1.0 or newer. There is no guarantee for backward or forward comparability. Please raise the issue on GitHub if something breaks.

The following R packages are required:

- RcppArmadillo (>= 14.2.2-1)
- Rcpp (>= 1.0.13-1)
- MASS (>= 7.3-60.2)
- irlba (>= 2.3.5.1)
- clue (>= 0.3-65)
- ciftiTools (>= 0.16.1)
- devtools

You can install them by running this code:

```r
if(!require(c("Rcpp","RcppArmadillo","MASS","irlba","clue","ciftiTools","devtools"))){
    install.packages(c("Rcpp","RcppArmadillo","MASS","irlba","clue","ciftiTools","devtools"))
}
```

Then you can install Sparse ICA from github with:

``` r
library(devtools)
install_github("thebrisklab/SparseICA")
# Load the package
library(SparseICA)
```

## Tutorial

The `sparseICA()` is the main function of our Sparse ICA algorithm. It is implemented in both pure R and Rcpp.

### Explanation of Arguments  
```
sparseICA(
    xData, n.comp, nu = "BIC", nu_list = seq(0.1, 4, 0.1), U.list = NULL,
    whiten = c('eigenvec', 'sqrtprec', 'none'), lngca = FALSE,
    orth.method = c('svd', 'givens'), method = c("C", "R"),
    restarts = 40, use_irlba = TRUE, eps = 1e-06, maxit = 500,
    verbose = TRUE, BIC_verbose = FALSE, converge_plot = FALSE, col.stand = TRUE,
    row.stand = FALSE, iter.stand = 5, positive_skewness = TRUE
)
```
- `xData`: A numeric matrix of input data with dimensions P x T, where P is the number of features and T is the number of samples.
- `n.comp`: An integer specifying the number of components to estimate.
- `nu`: A positive numeric value or a character "BIC" specifying the tuning parameter controlling the balance between accuracy and sparsity of the results. It can be selected using a BIC-like criterion (\code{"BIC"}) or based on expert knowledge (a positive number). Default is "BIC".
- `nu_list`: A numeric vector specifying the list of candidate tuning parameters. Default is \code{seq(0.1, 4, 0.1)}.
- `U.list`: An optional matrix specifying the initialization of the U matrix. Default is \code{NULL}.
- `whiten`: A character string specifying the method for whitening the input \code{xData}. Options are \code{"eigenvec"}, \code{"sqrtprec"}, \code{"lngca"}, or \code{"none"}. Default is \code{"eigenvec"}.
- `lngca`: A logical value indicating whether to perform Linear Non-Gaussian Component Analysis (LNGCA). Default is \code{FALSE}.
- `orth.method`: A character string specifying the method used for generating initial values for the U matrix. Default is \code{"svd"}.
- `method`: A character string specifying the computation method. If \code{"C"} (default), C code is used for most computations for better performance. If \code{"R"}, computations are performed entirely in R.
- `restarts`: An integer specifying the number of random initializations for optimization. Default is 40.
- `use_irlba`: A logical value indicating whether to use the \code{irlba} method for fast truncated Singular Value Decomposition (SVD) during whitening. This can improve memory efficiency for intermediate datasets. Default is \code{TRUE}.
- `eps`: A numeric value specifying the convergence threshold. Default is \code{1e-6}.
- `maxit`: An integer specifying the maximum number of iterations for the Sparse ICA method using Laplace density. Default is 500.
- `verbose`: A logical value indicating whether to print convergence information during execution. Default is \code{TRUE}.
- `BIC_verbose`: A logical value indicating whether to print BIC selection information. Default is \code{FALSE}.
- `converge_plot`: A logical value indicating whether to generate a line plot showing the convergence trace. Default is \code{FALSE}.
- `col.stand`: A logical value indicating whether to standardize columns. For each column, the mean of the entries in the column equals 0, and the variance of the entries in the column equals 1. Default is \code{TRUE}.
- `row.stand`: A logical value indicating whether to standardize rows. For each row, the mean of the entries in the row equals 0, and the variance of the entries in the row equals 1. Default is \code{FALSE}.
- `iter.stand`: An integer specifying the number of iterations for achieving both row and column standardization when \code{col.stand = TRUE} and \code{row.stand = TRUE}. Default is 5.
- `positive_skewness`: A logical value indicating whether to enforce positive skewness on the estimated components. Default is \code{TRUE}.

### Explanation of Output
The output will be a list with the following components as such:
- `loglik`: The minimal log-likelihood value among the random initializations.
- `estS`: A numeric matrix of estimated sparse independent components with dimensions P x Q.
- `estM`: The estimated mixing matrix with dimensions Q x T.
- `estU`: The estimated U matrix with dimensions Q x Q.
- `whitener`: The whitener matrix used for data whitening.
- `converge`: The convergence trace of the difference norm of the U matrix.
- `BIC`: A numeric vector of BIC values corresponding to each candidate `nu` in `nu_list`.
- `nu_list`: A numeric vector of candidate tuning parameter values.
- `best_nu`: The optimal nu selected based on the BIC-like criterion.


## Example
Load the example data:
``` r
data(example_sim123)
```

### 1. Visualization of example data
- The true "source signals" - `smat`:
```r
par(mfrow=c(1,3))
image(matrix(-example_sim123$smat[,1],33))
image(matrix(-example_sim123$smat[,2],33))
image(matrix(-example_sim123$smat[,3],33))
```
<img src="fig/true123.png" width="1200" />

- The true time series in the mixing matrix - `mmat`:
```r
par(mfrow=c(3,1))
plot(example_sim123$mmat[1,],type = "l",xlab = "Time",ylab = "")
plot(example_sim123$mmat[2,],type = "l",xlab = "Time",ylab = "")
plot(example_sim123$mmat[3,],type = "l",xlab = "Time",ylab = "")
```
<img src="fig/trueM.png" width="1200" />

- The simulated data matrix at time points T = 12, 23, and 35 - `xmat`:
```r
par(mfrow=c(1,3))
image(matrix(example_sim123$xmat[,12],33))
image(matrix(example_sim123$xmat[,23],33))
image(matrix(example_sim123$xmat[,35],33))
```
<img src="fig/xmat.png" width="1200" />


### 2. Run Sparse ICA
- Perform our Sparse ICA algorithm using `sparseICA` function. The tuning parameter `nu` is selected using the BIC-like criterion.
```r
my_sparseICA = sparseICA(xData = example_sim123$xmat, n.comp = 3, nu = "BIC", method = "C",restarts = 40,
                           eps = 1e-6, maxit = 500, verbose=TRUE)
```

<img src="fig/BIC.png" width="1200" />

The selected optimal `nu` is 1.2.

### 3. Visualization of Sparse ICA results
- Match estimates with the truth.
```r
matched_res=matchICA(my_sparseICA$estS,smat,example_sim123$smat)
```

- Check the estimated source signals.
```r
par(mfrow=c(1,3))
image(matrix(-matched_res$estS[,1],33,33))
image(matrix(-matched_res$estS[,2],33,33))
image(matrix(-matched_res$estS[,3],33,33))
```
<img src="fig/estS.png" width="1200" />

- Check the estimated time series in the mixing matrix.
```r
par(mfrow=c(3,1))
plot(matched_res$estM[1,],type = "l",xlab = "Time",ylab = "")
plot(matched_res$estM[2,],type = "l",xlab = "Time",ylab = "")
plot(matched_res$estM[3,],type = "l",xlab = "Time",ylab = "")
```
<img src="fig/estM.png" width="1200" />

- Check the correlations.
```r
> cor(smat[,1],matched_res$estS[,1])
[1] 0.9974975
> cor(smat[,2],matched_res$estS[,2])
[1] 0.9975088
> cor(smat[,3],matched_res$estS[,3])
[1] 0.9900951

> cor(mmat[1,],matched_res$estM[1,])
[1] 0.9645858
> cor(mmat[2,],matched_res$estM[2,])
[1] 0.99116
> cor(mmat[3,],matched_res$estM[3,])
[1] 0.9924354
```

## Citation
Those using the **SparseICA** software should cite:    
[Wang Z., Gaynanova, I., Aravkin, A., Risk, B. B. (2024). Sparse Independent Component Analysis with an Application to Cortical Surface fMRI Data in Autism. Journal of the American Statistical Association, 119(548), 2508-2520.](https://doi.org/10.1080/01621459.2024.2370593)

