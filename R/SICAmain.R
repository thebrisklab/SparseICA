
#' Function for generating random starting points
#' @param p The number of rows.
#' @param d The number of columns.
#' @param runs The number of random starts.
#' @param orth.method The method used for generating initial values of U matrix. The default is "svd".
#' @return A list of random initialization of matrices.
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method)
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]
    } else {
      temp = matrix(rnorm(p*d),p,d)
      W.list[[i]] <- svd(temp)$u
    }
  }
  return(W.list)
}


#' The function for perform whitening.
#'
#' @param X The data matrix with dimension P x T.
#' @param n.comp The number of components.
#' @param center.row Whether to center the row of data. Default is FALSE.
#' @param use_irlba Whether to use the irlba method to perform fast truncated singular value decomposition in whitening step, helpful for memorying intermediate dataset. Default is TRUE.
#' @return A list including the whitener matrix, the whitened data matrix, and the mean of the input data.
#'
#' @import MASS
#' @import irlba
#' 
#' @export
whitener <- function(X, n.comp = ncol(X), center.row = FALSE, use_irlba = TRUE) {
  if (ncol(X) > nrow(X)) warning("X is whitened with respect to columns")
  
  # Centering the data
  x.center <- scale(X, center = TRUE, scale = FALSE)
  if (center.row) x.center <- x.center - rowMeans(x.center)
  
  # Singular Value Decomposition (SVD) with `irlba`
  n.rep <- nrow(x.center)
  if(use_irlba) {
    svd.x <- irlba::irlba(x.center, nu = n.comp, nv = n.comp)
  } else {
    svd.x <- svd(x.center, nu = n.comp, nv = n.comp)
  }
  
  # Whitened data and parameters
  whitener_matrix <- t(ginv(svd.x$v %*% diag(svd.x$d[1:n.comp]) / sqrt(n.rep - 1)))
  Z <- sqrt(n.rep - 1) * svd.x$u
  mean_vec <- colMeans(X)
  
  list(whitener = whitener_matrix, Z = Z, mean = mean_vec)
}


#' Soft-threshold function
#' 
#' @param x The input scalar.
#' @param nu The tuning parameter.
#' @param lambda The lambda parameter of the Laplace density.
soft_thresh_R = function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin = pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}


#' Relax-and-split ICA Function for Sparse ICA wrapper
#'
#' This function performs Sparse Independent Component Analysis (Sparse ICA), implemented in both pure R and RCpp for efficiency.
#'
#' @param xData A numeric matrix of input data with dimensions P x T, where P is the number of features and T is the number of samples.
#' @param n.comp An integer specifying the number of components to estimate.
#' @param nu A numeric tuning parameter controlling the balance between accuracy and sparsity of the results. It can be selected using a BIC-like criterion or based on expert knowledge. Default is 1.
#' @param U.list An optional matrix specifying the initialization of the U matrix. Default is \code{NULL}.
#' @param whiten A character string specifying the method for whitening the input \code{xData}. Options are \code{"eigenvec"}, \code{"sqrtprec"}, \code{"lngca"}, or \code{"none"}. Default is \code{"eigenvec"}.
#' @param lngca A logical value indicating whether to perform Linear Non-Gaussian Component Analysis (LNGCA). Default is \code{FALSE}.
#' @param orth.method A character string specifying the method used for generating initial values for the U matrix. Default is \code{"svd"}.
#' @param method A character string specifying the computation method. If \code{"C"} (default), C code is used for most computations for better performance. If \code{"R"}, computations are performed entirely in R.
#' @param restarts An integer specifying the number of random initializations for optimization. Default is 40.
#' @param use_irlba A logical value indicating whether to use the \code{irlba} method for fast truncated Singular Value Decomposition (SVD) during whitening. This can improve memory efficiency for intermediate datasets. Default is \code{TRUE}.
#' @param eps A numeric value specifying the convergence threshold. Default is \code{1e-6}.
#' @param maxit An integer specifying the maximum number of iterations for the Sparse ICA method using Laplace density. Default is 500.
#' @param verbose A logical value indicating whether to print convergence information during execution. Default is \code{FALSE}.
#' @param converge_plot A logical value indicating whether to generate a line plot showing the convergence trace. Default is \code{FALSE}.
#' @param col.stand A logical value indicating whether to standardize columns. For each column, the mean of the entries in the column equals 0, and the variance of the entries in the column equals 1. Default is \code{TRUE}.
#' @param row.stand A logical value indicating whether to standardize rows. For each row, the mean of the entries in the row equals 0, and the variance of the entries in the row equals 1. Default is \code{FALSE}.
#' @param iter.stand An integer specifying the number of iterations for achieving both row and column standardization when \code{col.stand = TRUE} and \code{row.stand = TRUE}. Default is 5.
#' @param positive_skewness A logical value indicating whether to enforce positive skewness on the estimated components. Default is \code{TRUE}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{loglik}}{The minimal log-likelihood value among the random initializations.}
#'   \item{\code{estS}}{A numeric matrix of estimated sparse independent components with dimensions P x Q.}
#'   \item{\code{estU}}{The estimated U matrix with dimensions Q x Q.}
#'   \item{\code{estM}}{The estimated mixing matrix with dimensions Q x T.}
#'   \item{\code{whitener}}{The whitener matrix used for data whitening.}
#'   \item{\code{converge}}{Convergence information for the U matrix.}
#' }
#'
#' @import MASS
#' @import irlba
#' @import Rcpp
#' @import RcppArmadillo
#'
#' @export
relax_and_split_ICA <- function(
    xData, n.comp, nu = 1, U.list = NULL,
    whiten = c('eigenvec', 'sqrtprec', 'none'), lngca = FALSE,
    orth.method = c('svd', 'givens'), method = c("C", "R"),
    restarts = 40, use_irlba = TRUE, eps = 1e-06, maxit = 500,
    verbose = FALSE, converge_plot = FALSE, col.stand = TRUE,
    row.stand = FALSE, iter.stand = 5, positive_skewness = TRUE
) {
  #start.time <- Sys.time()
  
  # Argument parsing and initialization
  xData <- as.matrix(xData)
  whiten <- match.arg(whiten)
  orth.method <- match.arg(orth.method)
  method <- match.arg(method)
  
  if (verbose) message("Centering/scaling data.")
  
  xData_centered = scale(xData, center=TRUE, scale=FALSE)
  # Center and optionally scale data
  xData <- scale(xData, center = TRUE, scale = FALSE)
  if (col.stand || row.stand) {
    for (i in seq_len(iter.stand)) {
      if (col.stand) xData <- scale(xData)
      if (row.stand) xData <- t(scale(t(xData)))
    }
  }
  
  p <- ncol(xData)
  d <- n.comp
  t <- nrow(xData)
  
  if (d > p) stop("Number of components (n.comp) must not exceed number of columns in data.")
  
  # Whitening
  if (verbose) message("Whitening data.")
  whitener <- diag(d) # Default to identity
  
  if (lngca) {
    temp <- whitener(X = xData, n.comp = t - 1, use_irlba = use_irlba)
    xData <- temp$Z
    whitener <- temp$whitener
  } else if (whiten == "eigenvec") {
    temp <- whitener(X = xData, n.comp = d, use_irlba = use_irlba)
    xData <- temp$Z
    whitener <- temp$whitener
  } else if (whiten == "sqrtprec") {
    evd <- svd(cov(xData))
    whitener <- evd$u %*% diag(1 / sqrt(evd$d)) %*% t(evd$u)
    xData <- xData %*% whitener
  }
  
  # Initialize U.list if not provided
  if (is.null(U.list)) {
    U.list <- gen.inits(p = d, d = d, runs = restarts, orth.method = orth.method)
  }
  runs <- length(U.list)
  
  if (verbose) message("Starting Sparse ICA with ",runs," initialization(s).")
  
  ##############################################################################
  ############################ R implementation ################################
  ##############################################################################
  
  if (method=="R"){
    
    # Create a NULL list for storing outcomes
    out_list <- NULL
    
    converge = c()
    converge_plot = converge_plot
    
    # Store log likelihood
    loglik = c()

    for (k in 1:runs) { 
      # Initial value of V(0)
      newV = xData %*% U.list[[k]]
      lagU = diag(nrow = d)
      
      # loop for relax_laplace
      for (i in 1:maxit) {
        iteration = i
        # update U:
        txv = t(xData)%*%newV
        svd.txv = La.svd(txv)
        newU = svd.txv$u%*%svd.txv$vt
        # update V:
        newV = xData%*%newU
        newV = soft_thresh_R(newV, nu = nu, lambda = sqrt(2)/2)
        converge[iteration] = max(Mod(Mod(diag(newU %*% t(lagU))) - 1))
        lagU = newU
        if (converge[iteration] < eps) {break}
      }# End the loop for relax_laplace
      
      lambda = sqrt(2)/2
      out_list[[k]]=list()
      loglik[k] = sum(-abs(newV)/lambda-log(2*lambda)) - 1/(2*nu)*(norm(newV-xData%*%newU,type = "F"))^2
      out_list[[k]]$loglik = loglik[k]
      out_list[[k]]$estS = newV
      out_list[[k]]$estU = newU
      out_list[[k]]$converge = converge
      
    } #-------- end the loop to make a list
    
    final_out_list = out_list[[which.max(loglik)]]# update this
    final_out_list$whitener = whitener
    final_out_list$estM = est.M.ols(final_out_list$estS,xData_centered)
    if(positive_skewness){
      temp <- signchange(t(final_out_list$estS),t(final_out_list$estM))
      final_out_list$estS <- t(temp$S)
      final_out_list$estM <- t(temp$M)
    }
    
    if(verbose){
      for (converge_i in 1:length(final_out_list$converge)) {
        message("Iteration ",converge_i,", tol = ",final_out_list$converge[converge_i],".")}
      }
    
    if (converge_plot) {
      plot(final_out_list$converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
    }
  }else if (method=="C"){

    # Create a NULL list for storing outcomes
    final_out_list = NULL
    
    final_Rcpp = runs_relax_laplace(xData,U.list,runs,n.comp,nu,sqrt(2)/2,maxit,eps)
    
    if(verbose){
      converge = final_Rcpp$converge[final_Rcpp$converge!=0]
      for (converge_i in 1:length(converge)) {
        message("Iteration ",converge_i,", tol = ",converge[converge_i],".")
      }
    }

    final_out_list$loglik = final_Rcpp$loglik
    final_out_list$estS = final_Rcpp$newV
    final_out_list$estU = final_Rcpp$newU
    final_out_list$converge = final_Rcpp$converge
    final_out_list$whitener = whitener
    final_out_list$estM = est.M.ols(final_out_list$estS,xData_centered)
    if(positive_skewness){
      temp <- signchange(t(final_out_list$estS),t(final_out_list$estM))
      final_out_list$estS <- t(temp$S)
      final_out_list$estM <- t(temp$M)
    }

    if (converge_plot) {
      plot(final_out_list$converge[which(final_out_list$converge!=0)], main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
    }
    
  } # end of RCpp implementation
    return(final_out_list)
}


#' BIC-like Criterion for Tuning Parameter Selection in Sparse ICA
#'
#' This function uses a BIC-like criterion to select the optimal tuning parameter \code{nu} for Sparse ICA.
#'
#' @param xData A numeric matrix of input data with dimensions P x T, where P is the number of features and T is the number of samples.
#' @param n.comp An integer specifying the number of components to estimate.
#' @param nu_list A numeric vector specifying the list of candidate tuning parameters. Default is \code{seq(0.1, 4, 0.1)}.
#' @param whiten A character string specifying the method for whitening the input \code{xData}. Options are \code{"eigenvec"}, \code{"sqrtprec"}, or \code{"none"}. Default is \code{"eigenvec"}.
#' @param lngca A logical value indicating whether to perform Linear Non-Gaussian Component Analysis (LNGCA). Default is \code{FALSE}.
#' @param orth.method A character string specifying the method for generating initial values of the U matrix. Default is \code{"svd"}.
#' @param method A character string specifying the computation method. If \code{"C"} (default), C code is used for Sparse ICA to improve performance. If \code{"R"}, computations are performed entirely in R.
#' @param use_irlba A logical value indicating whether to use the \code{irlba} method for fast truncated Singular Value Decomposition (SVD) during whitening. This can improve memory efficiency for intermediate datasets. Default is \code{TRUE}.
#' @param eps A numeric value specifying the convergence threshold. Default is \code{1e-6}.
#' @param maxit An integer specifying the maximum number of iterations for the Sparse ICA method using Laplace density. Default is 500.
#' @param verbose A logical value indicating whether to print convergence information during execution. Default is \code{FALSE}.
#' @param col.stand A logical value indicating whether to standardize columns. For each column, the mean of the entries in the column equals 0, and the variance of the entries in the column equals 1. Default is \code{TRUE}.
#' @param row.stand A logical value indicating whether to standardize rows. For each row, the mean of the entries in the row equals 0, and the variance of the entries in the row equals 1. Default is \code{FALSE}.
#' @param iter.stand An integer specifying the number of iterations for achieving both row and column standardization when \code{col.stand = TRUE} and \code{row.stand = TRUE}. Default is 5.
#' @param BIC_plot A logical value indicating whether to generate a plot showing the trace of BIC values for different \code{nu} candidates. Default is \code{FALSE}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{BIC}}{A numeric vector of BIC values corresponding to each candidate \code{nu} in \code{nu_list}.}
#'   \item{\code{nu_list}}{A numeric vector of candidate tuning parameter values.}
#'   \item{\code{best_nu}}{The optimal \code{nu} selected based on the BIC-like criterion.}
#' }
#'
#' @import MASS
#' @import irlba
#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @examples
#' \donttest{
#' #get simulated data
#' data(example_sim123)
#'
#' select_sparseICA = BIC_sparseICA(xData = example_sim123$xmat, n.comp = 3, 
#'       method="C", BIC_plot = TRUE,verbose = TRUE, nu_list = seq(0.1,4,0.1))
#'
#' (my_nu = select_sparseICA$best_nu)
#' }
#' 
#' @export
BIC_sparseICA <- function(
    xData, n.comp, nu_list = seq(0.1, 4, 0.1),
    whiten = c('eigenvec', 'sqrtprec', 'none'), lngca = FALSE,
    orth.method = c('svd', 'givens'), method = c("C", "R"), 
    use_irlba = TRUE, eps = 1e-06, maxit = 500, verbose = FALSE,
    col.stand = TRUE, row.stand = FALSE, iter.stand = 0, BIC_plot = FALSE
) {
  start.time <- Sys.time()
  
  # Ensure input matrix format and parse arguments
  xData <- as.matrix(xData)
  whiten <- match.arg(whiten)
  orth.method <- match.arg(orth.method)
  method <- match.arg(method)
  
  if (verbose) message("Centering and scaling data.")
  
  # Center data
  xData_centered <- scale(xData, center = TRUE, scale = FALSE)
  v <- nrow(xData)
  t <- ncol(xData)
  
  if (verbose) message("Initializing reference Sparse ICA.")
  
  # Run reference Sparse ICA with very small `nu`
  ref_sparseICA <- relax_and_split_ICA(
    xData = xData, n.comp = n.comp, whiten = whiten,
    orth.method = orth.method, method = method, restarts = 1,
    use_irlba = use_irlba, nu = 1e-10, eps = eps, maxit = maxit,
    verbose = verbose, converge_plot = FALSE,
    col.stand = col.stand, row.stand = row.stand, iter.stand = iter.stand
  )
  
  my_W_list <- list(ref_sparseICA$estU) # Warm start for iterations
  
  # Initialize storage for BIC values
  out_BIC <- numeric(length(nu_list))
  
  if (verbose) message("Starting BIC calculation over ",length(nu_list)," nu values.")
  
  # Loop through `nu_list`
  for (i in seq_along(nu_list)) {
    my_nu <- nu_list[i]
    
    temp_sparseICA <- relax_and_split_ICA(
      xData = xData, n.comp = n.comp, U.list = my_W_list,
      whiten = whiten, orth.method = orth.method, method = method,
      restarts = 1, nu = my_nu, eps = eps, maxit = maxit,
      verbose = verbose, converge_plot = FALSE, use_irlba = use_irlba,
      col.stand = col.stand, row.stand = row.stand, iter.stand = iter.stand
    )
    
    S_nu <- temp_sparseICA$estS
    ginv_S <- ginv(crossprod(S_nu))
    residuals <- xData_centered - S_nu %*% (ginv_S %*% t(S_nu) %*% xData_centered)
    e_nu <- sum(residuals^2)
    
    # Calculate BIC
    sparsity_term <- sum(S_nu != 0) * log(v * t) / (v * t)
    out_BIC[i] <- log(e_nu / (v * t)) + sparsity_term
    
    # Update warm start
    my_W_list <- list(temp_sparseICA$estU)
  }
  
  # Select best `nu` based on minimum BIC
  best_nu <- nu_list[which.min(out_BIC)]
  
  if (verbose) message("The best nu selected by BIC is ", best_nu, ".")
  
  # Optional BIC plot
  if (BIC_plot) {
    plot(nu_list,out_BIC, main = "BIC Plot", xlab = "nu", ylab = "BIC", type = "l")
    abline(v=best_nu,col="red")
  }
  
  # Compile results
  result <- list(
    BIC = out_BIC,
    nu_list = nu_list,
    best_nu = best_nu
  )
  
  if (verbose) message("BIC selection completed in ",round(difftime(Sys.time(),start.time,units = "secs"),2)," seconds.")
  return(result)
}

#' Sparse Independent Component Analysis (Sparse ICA) Function
#'
#' This function performs Sparse Independent Component Analysis (Sparse ICA), implemented in both pure R and RCpp for efficiency.
#'
#' @param xData A numeric matrix of input data with dimensions P x T, where P is the number of features and T is the number of samples.
#' @param n.comp An integer specifying the number of components to estimate.
#' @param nu A positive numeric value or a character "BIC" specifying the tuning parameter controlling the balance between accuracy and sparsity of the results. It can be selected using a BIC-like criterion (\code{"BIC"}) or based on expert knowledge (a positive number). Default is "BIC".
#' @param nu_list A numeric vector specifying the list of candidate tuning parameters. Default is \code{seq(0.1, 4, 0.1)}.
#' @param U.list An optional matrix specifying the initialization of the U matrix. Default is \code{NULL}.
#' @param whiten A character string specifying the method for whitening the input \code{xData}. Options are \code{"eigenvec"}, \code{"sqrtprec"}, \code{"lngca"}, or \code{"none"}. Default is \code{"eigenvec"}.
#' @param lngca A logical value indicating whether to perform Linear Non-Gaussian Component Analysis (LNGCA). Default is \code{FALSE}.
#' @param orth.method A character string specifying the method used for generating initial values for the U matrix. Default is \code{"svd"}.
#' @param method A character string specifying the computation method. If \code{"C"} (default), C code is used for most computations for better performance. If \code{"R"}, computations are performed entirely in R.
#' @param restarts An integer specifying the number of random initializations for optimization. Default is 40.
#' @param use_irlba A logical value indicating whether to use the \code{irlba} method for fast truncated Singular Value Decomposition (SVD) during whitening. This can improve memory efficiency for intermediate datasets. Default is \code{TRUE}.
#' @param eps A numeric value specifying the convergence threshold. Default is \code{1e-6}.
#' @param maxit An integer specifying the maximum number of iterations for the Sparse ICA method using Laplace density. Default is 500.
#' @param verbose A logical value indicating whether to print convergence information during execution. Default is \code{TRUE}.
#' @param BIC_verbose A logical value indicating whether to print BIC selection information. Default is \code{FALSE}.
#' @param converge_plot A logical value indicating whether to generate a line plot showing the convergence trace. Default is \code{FALSE}.
#' @param col.stand A logical value indicating whether to standardize columns. For each column, the mean of the entries in the column equals 0, and the variance of the entries in the column equals 1. Default is \code{TRUE}.
#' @param row.stand A logical value indicating whether to standardize rows. For each row, the mean of the entries in the row equals 0, and the variance of the entries in the row equals 1. Default is \code{FALSE}.
#' @param iter.stand An integer specifying the number of iterations for achieving both row and column standardization when \code{col.stand = TRUE} and \code{row.stand = TRUE}. Default is 5.
#' @param positive_skewness A logical value indicating whether to enforce positive skewness on the estimated components. Default is \code{TRUE}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{loglik}}{The minimal log-likelihood value among the random initializations.}
#'   \item{\code{estS}}{A numeric matrix of estimated sparse independent components with dimensions P x Q.}
#'   \item{\code{estM}}{The estimated mixing matrix with dimensions Q x T.}
#'   \item{\code{estU}}{The estimated U matrix with dimensions Q x Q.}
#'   \item{\code{whitener}}{The whitener matrix used for data whitening.}
#'   \item{\code{converge}}{The trace of convergence for the U matrix.}
#'   \item{\code{BIC}}{A numeric vector of BIC values corresponding to each candidate \code{nu} in \code{nu_list}.}
#'   \item{\code{nu_list}}{A numeric vector of candidate tuning parameter values.}
#'   \item{\code{best_nu}}{The optimal \code{nu} selected based on the BIC-like criterion.}
#' }
#'
#' @import MASS
#' @import irlba
#' @import Rcpp
#' @import RcppArmadillo
#' @examples
#' \donttest{
#' #get simulated data
#' data(example_sim123)
#'
#' my_sparseICA <- sparseICA(xData = example_sim123$xmat, n.comp = 3, nu = "BIC", method="C",
#'       restarts = 40, eps = 1e-6, maxit = 500, verbose=TRUE)
#'
#' res_matched <- matchICA(my_sparseICA$estS,example_sim123$smat)
#'
#' # Visualize the estimated components
#' oldpar <- par()$mfrow
#' par(mfrow=c(1,3))
#' image(matrix(res_matched[,1],33,33))
#' image(matrix(res_matched[,2],33,33))
#' image(matrix(res_matched[,3],33,33))
#' par(mfrow=oldpar)
#' }
#' 
#' @export
sparseICA <- function(
    xData, n.comp, nu = "BIC", nu_list = seq(0.1, 4, 0.1), U.list = NULL,
    whiten = c('eigenvec', 'sqrtprec', 'none'), lngca = FALSE,
    orth.method = c('svd', 'givens'), method = c("C", "R"),
    restarts = 40, use_irlba = TRUE, eps = 1e-06, maxit = 500,
    verbose = TRUE, BIC_verbose = FALSE, converge_plot = FALSE, col.stand = TRUE,
    row.stand = FALSE, iter.stand = 5, positive_skewness = TRUE
) {
  start.time1 <- Sys.time()
  
  if(nu == "BIC"){
    if (verbose) message("Selecting the optimal nu using BIC-like criterion.")
    BIC_result <- BIC_sparseICA(
      xData = xData, n.comp = n.comp, nu_list = nu_list,
      whiten = whiten, lngca = lngca, orth.method = orth.method,
      method = method, use_irlba = use_irlba, eps = eps, maxit = maxit,
      verbose = BIC_verbose, col.stand = col.stand, row.stand = row.stand,
      iter.stand = iter.stand, BIC_plot = FALSE)
    nu <- BIC_result$best_nu
    if (verbose) message("The optimal nu is ",nu,".")
    
    sparseICA_result <- relax_and_split_ICA(
      xData = xData, n.comp = n.comp, U.list = U.list,
      whiten = whiten, lngca = lngca, orth.method = orth.method,
      method = method, restarts = restarts, use_irlba = use_irlba,
      nu = nu, eps = eps, maxit = maxit, verbose = verbose,
      converge_plot = converge_plot, col.stand = col.stand,
      row.stand = row.stand, iter.stand = iter.stand,
      positive_skewness = positive_skewness
    )
    
    sparseICA_result$best_nu <- BIC_result$best_nu
    sparseICA_result$BIC <- BIC_result$BIC
    sparseICA_result$nu_list <- BIC_result$nu_list
    
    if (verbose) message("SparseICA is completed in ",round(difftime(Sys.time(),start.time1,units = "secs"),2)," seconds.")
    
    return(sparseICA_result)
  } else {
    sparseICA_result <- relax_and_split_ICA(
      xData = xData, n.comp = n.comp, U.list = U.list,
      whiten = whiten, lngca = lngca, orth.method = orth.method,
      method = method, restarts = restarts, use_irlba = use_irlba,
      nu = nu, eps = eps, maxit = maxit, verbose = verbose,
      converge_plot = converge_plot, col.stand = col.stand,
      row.stand = row.stand, iter.stand = iter.stand,
      positive_skewness = positive_skewness
    )
    
    if (verbose) message("SparseICA is completed in ",round(difftime(Sys.time(),start.time1,units = "secs"),2)," seconds.")
    
    return(sparseICA_result)
  }
}

