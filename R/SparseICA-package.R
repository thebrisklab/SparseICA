## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline
#' @importFrom stats coef cov lm rnorm runif
#' @useDynLib SparseICA, .registration = TRUE
## usethis namespace: end
NULL


#' Sign change for S matrix to image
#' @param S S matrix, P x Q.
#' @param M M matrix, Q x T.
#'
#' @return  a list of positive S and positive M.
#' @export
signchange = function(S,M=NULL) {
  signskew = sign(apply(S,1,function(x) mean(x^3)))
  newS = signskew*S        # t(signskew*t(S))
  if(!is.null(M)) {
    newM = t(signskew*t(M))
  }else {newM=M}
  return(list(S=newS,M=newM))
}



#' Function for generating random starting points
#' @param p Number of rows
#' @param d Number of columns
#' @param runs the number of random starts
#' @param orth.method orthodox method
#' @return a list of initialization of matrices.
#' @examples gen.inits(2,3,3,"svd")
#' @export
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
  W.list
}

#' For a given angle theta, returns a d x d Givens rotation matrix
#' @param theta the value of theta.
#' @param d the value of d.
#' @param which the value of which.
givens.rotation <- function(theta=0, d=2, which=c(1,2))
{
  # For a given angle theta, returns a d x d Givens rotation matrix
  #
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}


#' Convert angle vector into orthodox matrix
#'
#' @param theta vector of angles theta
#' @return an orthodox matrix
theta2W = function(theta){
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W
      index = index + 1
    }
  }
  W
}



#' Whitening Function
#'
#' @param X The data matrix P x T
#' @param n.comp The number of components
#' @param center.row Whether center the row of data
#' @param irlba Whether use the irlba package
#' @return a whitener matrix
#'
#' @import MASS
#' @import irlba
#' @export
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}


#' Estimate mixing matrix from estimates of components
#' @param sData S Dimension: P x Q
#' @param xData X Dimension: P x T
#' @param intercept default = TRUE
#' @return a mixing matrix M, dimension Q x T.
#' @export
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}


#' Soft-thresholding to update V
#' @param x input scalar
#' @param nu the tuning parameter
#' @param lambda the lambda parameter of the Laplace density
soft_thresh_R = function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin = pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}


#' Sparse ICA function
#'
#' The function to perform Sparse ICA, which is implemented in both pure R and RCpp.
#'
#' @param xData Input data matrix with dimension P x T. P is the number of features. T is the number of samples.
#' @param n.comp The number of components.
#' @param nu the tuning parameter controlling the accuracy and sparsity of the results. Should be selected by the BIC-like criterion "BIC_sparseICA_R()" or expert knowledge. The default is 1.
#' @param U.list The initialization of U matrix. Default is "NULL".
#' @param whiten The method for whitening input xData. Could take "eigenvec", "sqrtprec", "lngca", and "none". The default is "eigenvec".
#' @param orth.method The method used for generating initial values of U matrix. The default is "svd".
#' @param method If method == "C" (default) then C code is used to perform most of the computations, which makes the algorithm run faster. If method == "R" then computations are done exclusively in R.
#' @param restarts The number of initial points. Default is 40.
#' @param lambda The scale parameter in Laplace density. The default is sqrt(2)/2 to make the default situation with unit variance.
#' @param irlba Whether use the irlba method to perform fast truncated singular value decomposition in whitening step. The default is FALSE.
#' @param eps The convergence threshold. The default is 1e-6.
#' @param maxit The maximum number of iterations in the Sparse ICA method with Laplace density. The default number is 500.
#' @param verbose Whether print the information about convergence. The default is FALSE.
#' @param converge_plot Whether make a convergence plot for Sparse ICA method. The default is FALSE.
#' @param col.stand Whether standardize the data matrix column-wise. Default if TRUE.
#' @param row.stand Whether standardize the data matrix row-wise. Default if FALSE.
#' @param iter.stand The number of standardization. Default is 0.
#' @param positive_skewness Whether to make the estimated components have positive skewness. Default is TRUE.
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{loglik_restarts}}{The value of minimal log-likelihood among random initializations.}
#'       \item{\code{estS}}{The estimated sparse independent components matrix with dimension P x Q.}
#'       \item{\code{estU}}{The estimated U matrix with dimension Q x Q.}
#'       \item{\code{converge}}{The convergence of the U matrix.}
#'       \item{\code{distribution}}{The density used in ICA method.}
#'       \item{\code{whitener}}{The whitener matrix used to perform data whitening.}
#'       \item{\code{estM}}{Estimated mixing matrix with dimension Q x T.}}
#'
#' @import MASS
#' @import irlba
#' @import Rcpp
#' @import RcppArmadillo
#' @export
#'
#' @examples
#' \donttest{
#' #get simulated data
#' data(example_sim123)
#'
#' my_sparseICA = sparseICA(xData = xmat, n.comp = 3, nu = 1.2, U.list=NULL,
#'     whiten = "eigenvec", orth.method = "svd", method="C",restarts = 40,
#'     lambda = sqrt(2)/2, eps = 1e-6, maxit = 500, verbose=TRUE)
#'
#' a=matchICA(my_sparseICA$estS,smat)
#'
#' par(mfrow=c(1,3))
#' image(matrix(a[,1],33,33))
#' image(matrix(a[,2],33,33))
#' image(matrix(a[,3],33,33))
#' par(mfrow=c(1,1))
#' }
sparseICA = function(xData,n.comp,nu = 1,U.list=NULL,whiten = c('eigenvec','sqrtprec','lngca','none'),orth.method=c('svd','givens'), method = c("C","R"),restarts = 40, lambda = sqrt(2)/2, irlba = FALSE, eps = 1e-06, maxit = 500, verbose=FALSE, converge_plot = FALSE,col.stand=TRUE, row.stand=FALSE, iter.stand=0,positive_skewness=TRUE){

  start.time = Sys.time()

  ##################################################################################### 0
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  method = match.arg(method)

  if (verbose){
    cat("Centering.\n")
  }
  xData_centered = scale(xData, center=TRUE, scale=FALSE)
  # optional iterative center+scale included (recommend 5 times)
  if (iter.stand==0){
    # center the input data
    xData = scale(xData, center=TRUE, scale=FALSE)
  }else{
    for (i in 1:iter.stand) {
      if (row.stand){
        # standardization at each voxel
        xData[which(xData[,1]!=0),] = t(scale(t(xData[which(xData[,1]!=0),])))
      }
      if (col.stand){
        # standardization at each time point
        xData = scale(xData)
      }
    }
  }

  #if(irlba) require(irlba)
  p = ncol(xData) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. p = number of column
  d = n.comp # Extract the components need for V(0)
  t = nrow(xData)

  if (verbose){
    cat("Whitening.\n")
  }
  # Data whiten
  if (d > p) stop('d must be less than or equal to p')
  if (whiten=='eigenvec') {
    # Use whitener=='eigenvec' so that restarts.dbyd initiates from the span of the first d eigenvectors.
    temp = whitener(X = xData,n.comp = d,irlba=irlba)
    xData = temp$Z
    whitener = temp$whitener
    rm(temp)
  }else if (whiten=='sqrtprec') {
    est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
    evd.sigma = svd(est.sigma)
    whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
    xData = xData%*%whitener
  }else if (whiten=='lngca'){
    temp = whitener(X = xData,n.comp = (t-1),irlba=irlba)
    xData = temp$Z
    whitener = temp$whitener
    rm(temp)
  }else {
    whitener = diag(d)
  }

  ##############################################################################
  ############################ R implementation ################################
  ##############################################################################

  if (method=="R"){

  # Randomly make an input for U(0), by default here used orth.method = "svd"
  if (is.null(U.list)){
    U.list = gen.inits(p=d, d=d, runs = restarts, orth.method=orth.method)
  }
  runs = length(U.list)

  # Create a NULL list for storing outcomes
  out.list = NULL

  converge = c()
  converge_plot = converge_plot

  # Store log likelihood
  loglik = c()

  ##################################################################################
  #m=0
  if (verbose){
    cat("Start Sparse ICA, using ",runs," different random initializations. The convergence criterion is ",eps,".\n")
  }
  for (k in 1:runs) {  # make a list if restarts > 1, we should update the index of U.list if restarts > 1
    # Initial value of V(0)
    newV = xData %*% U.list[[k]]
    #lagV = newV
    lagU = diag(nrow = d)

    # loop for relax_laplace
    for (i in 1:maxit) {
      iteration = i
      # update U:
      txv = t(xData)%*%newV
      svd.txv = La.svd(txv)
      newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U
      # update V:
      newV = xData%*%newU
      newV = soft_thresh_R(newV, nu = nu, lambda = lambda)
      # if end after thresholding, ensures that the solution has exact zeros.
      #converge[iteration] = norm(lagV-newV,"F")
      #converge[iteration] = mean((lagV-newV)^2)#/ncol(lagV) # Mean squared norm
      converge[iteration] = max(Mod(Mod(diag(newU %*% t(lagU))) - 1))
      #converge[iteration] = mean((lagU-newU)^2) # Mean squared norm
      #lagV = newV
      lagU = newU
      if (converge[iteration] < eps) {break}
    }
    # End the loop for relax_laplace

    #m=m+1
    # # Store Warning
    # if(iteration==maxit){
    #   #cat("Does not converage!\n")
    #   m=m+1
    # }

    # Store the results
    #if (k != 1) {out.list[[k]] = out.list[[1]]} # update this
    out.list[[k]]=list()
    loglik[k] = sum(-abs(newV)/lambda-log(2*lambda)) - 1/(2*nu)*(norm(newV-xData%*%newU,type = "F"))^2
    out.list[[k]]$loglik = loglik[k]
    out.list[[k]]$estS = newV
    out.list[[k]]$estU = newU
    #out.list[[k]]$xData_whitened = xData
    out.list[[k]]$converge = converge

    if(verbose){
      for (converge_i in 1:length(converge)) {
        cat("Iteration ",converge_i,", tol=",converge[converge_i],".\n")
      }
    }

  } #-------- end the loop to make a list

  final.out.list = out.list[[which.max(loglik)]]# update this
  final.out.list$distribution = "Laplace"
  final.out.list$whitener = whitener
  final.out.list$loglik_restarts = loglik
  final.out.list$estM = est.M.ols(final.out.list$estS,xData_centered)
  if(positive_skewness){
    temp <- signchange(t(final.out.list$estS),t(final.out.list$estM))
    final.out.list$estS <- t(temp$S)
    final.out.list$estM <- t(temp$M)
  }
  #out.list$estM = solve(out.list$estW) %*% ginv(out.list$whitener)

  if (converge_plot) {
    plot(final.out.list$converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
  }
  end.time = Sys.time()
  # Print running time
  #time.taken = proc.time() - start.time
  if (verbose){
    print(end.time - start.time)
  }
  #time.taken=end.time - start.time
  #cat("Running time:",time.taken,"s.\n")

  return(final.out.list)
  } # end of R implementation

  ##############################################################################
  ############################ RCpp implementation #############################
  ##############################################################################

  if (method=="C"){

    # Randomly make an input for U(0), by default here used orth.method = "svd"
    if (is.null(U.list)){
      U.list = gen.inits(p=d, d=d, runs = restarts, orth.method=orth.method)
    }
    runs = length(U.list)

    if (verbose){
      cat("Start Sparse ICA, using ",runs," different random initializations. The convergence criterion is ",eps,".\n")
    }
    ###########################################################################################
    # Create a NULL list for storing outcomes
    out.list = NULL

    #converge = c()
    #converge_plot = converge_plot

    # Store log likelihood
    #loglik = c()

    final_Rcpp=runs_relax_laplace(xData,U.list,runs,n.comp,nu,lambda,maxit,eps)

    # if(verbose){
    #   cat("MESSAGE: Sparse ICA converges (<",eps,") in",sum(final_Rcpp$converge!=0),"iterations using",runs,"different random initiations!\n")
    # }
    if(verbose){
      converge = final_Rcpp$converge[final_Rcpp$converge!=0]
      for (converge_i in 1:length(converge)) {
        cat("Iteration ",converge_i,", tol=",converge[converge_i],".\n")
      }
    }

    out.list$loglik = final_Rcpp$loglik
    out.list$estS = final_Rcpp$newV
    out.list$estU = final_Rcpp$newU
    out.list$converge = final_Rcpp$converge
    out.list$distribution = "Laplace"
    out.list$whitener = whitener
    out.list$estM = est.M.ols(out.list$estS,xData_centered)
    if(positive_skewness){
      temp <- signchange(t(out.list$estS),t(out.list$estM))
      out.list$estS <- t(temp$S)
      out.list$estM <- t(temp$M)
    }

    if (converge_plot) {
      plot(out.list$converge[which(out.list$converge!=0)], main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
    }

    # Print running time
    end.time = Sys.time()
    if (verbose){
      print(end.time - start.time)
    }
    #cat("Running time:",time.taken,"s.\n")

    return(out.list)

  } # end of RCpp implementation

} # end of sparseICA function


#' BIC-like selection criterion function
#'
#' Use BIC-like criterion to select tuning parameter nu in Sparse ICA.
#'
#' @param xData Input data matrix with dimension P x T. P is the number of features. t is the number of samples.
#' @param n.comp The number of components.
#' @param nu_list the list of candidate tuning parameter. Default is seq(0.1,4,0.1).
#' @param U.list The initialization of U matrix. Default is "NULL".
#' @param whiten The method for whitening input xData. Could take "eigenvec", "sqrtprec", "lngca", and "none". The default is "eigenvec".
#' @param orth.method The method used for generating initial values of U matrix. The default is "svd".
#' @param method If method == "C" (default) then C code is used for Sparse ICA. If method == "R" then computations are done exclusively in R.
#' @param lambda The scale parameter in Laplace density. The default is sqrt(2)/2 to make the default situation with unit variance.
#' @param irlba Whether use the irlba method to perform fast truncated singular value decomposition in whitening step. The default is FALSE.
#' @param eps The convergence threshold. The default is 1e-6.
#' @param maxit The maximum number of iterations in the Sparse ICA method with Laplace density. The default number is 500.
#' @param verbose Whether print the information about convergence. The default is FALSE.
#' @param col.stand Whether standardize the data matrix column-wise. Default if TRUE.
#' @param row.stand Whether standardize the data matrix row-wise. Default if FALSE.
#' @param iter.stand The number of standardization. Default is 0.
#' @param BIC_plot Whether make a convergence plot for BIC selection. The default is FALSE.
#'
#' @return Function outputs a list including the following:
#' \describe{
#'       \item{\code{BIC}}{The list of BIC values corresponding to each candidate in the nu_list.}
#'       \item{\code{nu_list}}{The list of candidate tuning parameter nu.}
#'       \item{\code{best_nu}}{The best nu selected by BIC-like criterion.}}
#'
#' @import MASS
#' @import irlba
#' @import Rcpp
#' @import RcppArmadillo
#' @export
#'
#' @examples
#' \donttest{
#' #get simulated data
#' data(example_sim123)
#'
#' select_sparseICA = BIC_sparseICA(xData = xmat, n.comp = 3,U.list=NULL,
#'     whiten = "eigenvec", orth.method="svd", method="C", lambda = sqrt(2)/2,
#'     eps = 1e-6,maxit = 500, BIC_plot = TRUE,verbose=FALSE,
#'     nu_list = seq(0.1,4,0.1))
#'
#' my_nu = select_sparseICA$best_nu
#' }
BIC_sparseICA = function(xData,n.comp,nu_list = seq(0.1,4,0.1),U.list=NULL,whiten = c('eigenvec','sqrtprec','lngca','none'), orth.method=c('svd','givens'), method = c("C","R"), lambda = sqrt(2)/2, irlba = FALSE, eps = 1e-06, maxit = 500, verbose=FALSE,col.stand=TRUE, row.stand=FALSE, iter.stand=0, BIC_plot = FALSE){
  start.time = Sys.time()

  ##############################################################################
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)

  ##########################################################################################

  # reference Sparse ICA
  best_nu = 0.0000000001
  ref_sparseICA = sparseICA(xData = xData, n.comp = n.comp,whiten = whiten,orth.method=orth.method,
                              method=method,restarts = 1,irlba = irlba,
                              lambda = lambda, nu = best_nu,eps = eps,
                              maxit = maxit, verbose = verbose,
                              converge_plot = F,col.stand=col.stand, row.stand=row.stand, iter.stand=iter.stand)

  # center raw data matrix
  xData_centered = scale(xData,center = T,scale = F)

  v = dim(xData)[1]   # the number of voxels
  t = dim(xData)[2]   # the number of time points
  my_W_list = list(ref_sparseICA$estU)   # initial U for next start

  ##############################################################################
  # start BIC process
  out_BIC = rep(NA,length(nu_list))
  out_ratio = rep(NA,length(nu_list))
  out_sparsity = rep(NA,length(nu_list))

  for (i in 1:length(nu_list)) {
    my_nu = nu_list[i]
    temp_sparseICA = sparseICA(xData = xData,n.comp = n.comp,U.list = my_W_list,
                                 whiten = whiten,orth.method=orth.method,method=method,
                                 restarts = 1, lambda = lambda,
                                 nu = my_nu,eps = eps,maxit = maxit,
                               verbose = verbose, converge_plot = F,
                                 col.stand=col.stand, row.stand=row.stand, iter.stand=iter.stand)
    S_nu = temp_sparseICA$estS
    # generalized inverse?
    temp_ginv = ginv(crossprod(S_nu))
    temp_mat = temp_ginv%*%t(S_nu)%*%xData_centered
    #temp_mat = solve(crossprod(S_nu), t(S_nu))%*%xmat_centered
    e_nu = sum((xData_centered-S_nu%*%temp_mat)^2) # Calculate e_nu

    my_W_list = list(temp_sparseICA$estU)   # Warm start for next iteration

    out_ratio[i] = log(e_nu/(v*t))
    out_sparsity[i] = sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
    out_BIC[i] = log(e_nu/(v*t))+sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
  }

  best_nu = nu_list[which(out_BIC==min(out_BIC))][1]
  if(verbose){
    cat("The best nu selected by BIC is ",best_nu,".\n")
  }

  if (BIC_plot) {
    plot(nu_list,out_BIC, main = "BIC Plot", xlab = "nu", ylab = "BIC", type = "l")
    abline(v=best_nu,col="red")
  }

  out.list=NULL
  #out.list$ratio = out_ratio
  #out.list$sparsity = out_sparsity
  out.list$BIC = out_BIC
  out.list$nu_list = nu_list
  out.list$best_nu = best_nu

  # Print running time
  end.time=Sys.time()
  if (verbose){
    print(end.time - start.time)
  }

  return(out.list)
}



#' Match ICA results based on L2 distances and Hungarian
#'
#' @param S loading variable matrix
#' @param template template for match
#' @param M subject score matrix
#'
#' @return the match result
#' @import clue
#' @export
matchICA=function(S,template,M=NULL) {
  #require(clue)
  n.comp=ncol(S)
  n.obs=nrow(S)
  if(n.comp>n.obs) warning('Input should be n x d')
  if(!all(dim(template)==dim(S))) warning('Template should be n x d')
  S = t(S)
  template = t(template)
  l2.mat1=matrix(NA,nrow=n.comp,ncol=n.comp)
  l2.mat2=l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i,j]=sum((template[i,]-S[j,])^2)/n.obs
      l2.mat2[i,j]=sum((template[i,]+S[j,])^2)/n.obs
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)

  s.perm=t(perm)%*%S
  if(!is.null(M)) {
    M.perm=t(M)%*%perm
    return(list(S=t(s.perm),M=t(M.perm)))
  }  else {
    t(s.perm)
  }
}


