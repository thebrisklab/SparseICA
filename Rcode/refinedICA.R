##########################################################################
# soft_thresh to update matrix in relax and split method
# Laplace distribution is used, the default nu is 1
soft_thresh <- function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin <- pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}


##########################################################################
# Updated 2020.04.27 Did not change the code, but the original objective function is equivalent to this one under pre-whitened data
# Updated 2020.04.01 to change back the equation (the original objective function should be correct)
# Updated 2020.03.31 to change the equation
# Newton's method to find minimum/maximum (root for the first derivative)
prox_newton <- function(v, nu = 1, itm = 50, tol= 1e-6) {
  
  # based on the function of logistic: G(y) = ln(1+exp(a*y))
  
  # initialization
  z = v
  a = -pi/sqrt(3.0)
  n = 0
  e = exp(a*v)
  g = z - v + 2*a*nu*e/(1.0 + e) - a*nu # it is the function that we need to find its root
  h = 1.0 + 4*a^2*nu*e/(1.0 + e)^2 # it is the first derivative of g
  # Additional note: the logistic density problem in Peng's ppt is to find the minimum of nu*log(1.0 + e) + 0.5*(z - v)^2
  
  while(abs(sum(g)) >= tol) {
    z = z - g/h # the Newton's Method equation
    
    # update gradient hessian
    n = n + 1;
    e = exp(a*z);
    g = z - v + 2*a*nu*e/(1.0 + e) - a*nu
    h = 1.0 + 4*a^2*nu*e/(1.0 + e)^2
    
    if(n >= itm) {
      message(paste0("Already achieved maximum iteration=",itm))
      break} # force to quit when more than the maximum iteration
    
  } # end of while loop
  
  return(z)
  
} # end of prox_newton function


##########################################################################
# Updated 2020.04.28 for (1) rename relax-softmax to relax-laplace; 
#                        (2) add the log likelihood value for all the methods; 
#                        (3) allow data whitening with choices of whiten = c('eigenvec','sqrtprec','none'); 
#                        (4) allow different initialization methods for the start up unmixing (W) matrix, orth.method=c('svd','givens');
#                        (5) allow restarts > 1 for relax-laplace and relax-logistic; (FastICAs are originally allowed by mlcaFP);
#                        (6) allow to plots the convergence trend for relax-laplace and relax-logistic, by default it is hidden; 
#                        (7) note: the output format is different from the previous one. 
# Updated 2020.04.06 for convergence rule
# Updated 2020.04.01 for eps
# Updated 2020.03.30 for whiten
refinedICA = function(xData = xmat, n.comp = 3, W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 1, method = "relax_laplace", lambda = sqrt(2)/2, nu = 1, irlba = FALSE, eps = 1e-06, converge_plot = F){
  # Define a function to use refined ICA method. 
  # method = c("relax_logistic", "relax_laplace", "FastICA_logistic", "FastICA_tanh")
  
  ########################### 0
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  if(irlba) require(irlba)
  p = ncol(xData) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. 
  d = n.comp # Extract the components need for V(0)
  
  # Data whiten
  xData <- scale(xData, center=TRUE, scale=FALSE)
  if (d > p) stop('d must be less than or equal to p')
  if (whiten=='eigenvec') {
    # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
    # span of the first d eigenvectors.
    temp = whitener(X = xData,n.comp = p,irlba=irlba)
    xData = temp$Z
    whitener = temp$whitener
    rm(temp)
  }  else if (whiten=='sqrtprec') {
    est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
    evd.sigma = svd(est.sigma)
    whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
    xData = xData%*%whitener
  }
  else {
    whitener = diag(p)
  }

  W.list = gen.inits(p=p, d=d, runs = restarts.pbyd, orth.method=orth.method) # Randomly make an input for V(0), by default here used orth.method = "svd"
  runs = length(W.list)
  out.list = NULL
  
  converge = c()
  converge_plot = converge_plot
  loglik = c()

  ########################### 1
  if (method == "relax_laplace") {
    
    for (k in 1:runs) { #---------------------------------------- make a list if restarts.pbyd > 1
      newV = xData %*% W.list[[k]] # To be more accurate, we should update the index of W.list if restarts.pbyd > 1 
      lagV = newV
      for (i in 1:500) { # ---------- loop for relax_laplace
        iteration = i
        # update U:
        txv = t(xData)%*%newV
        svd.txv = La.svd(txv)
        newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U
        # update V:
        newV = xData%*%newU
        newV = soft_thresh(newV, nu = nu, lambda = lambda)
        # if end after thresholding, ensures that the solution has exact zeros. 
        converge[iteration] = norm(lagV-newV) / ncol(lagV) # frobICA(S1=newV, S2=lagV) # BRISK note: normalize by number of columns
        lagV = newV
        if (converge[iteration] < eps) {break}
      } #--------------------------- end the loop for relax_laplace
      
      if (k != 1) {out.list[[k]] = out.list[[1]]}
      out.list[[k]]$loglik = sum(-abs(lagV))
      loglik[k] = out.list[[k]]$loglik
      out.list[[k]]$estS = newV
      out.list[[k]]$converge = converge
      
    } #----------------------------------------------------------- end the loop to make a list
    
    out.list = out.list[[which.max(loglik)]]
    out.list$distribution = "Laplace"
    out.list$whitener = whitener
    
    if (converge_plot == F) {
      out.list = out.list[-3]
    } else {
      converge = out.list$converge
      plot = recordPlot()
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
      out.list$converge = plot
    }
    
  }
  
  ########################### 2
  if (method == "relax_logistic") {
    
    for (k in 1:runs) { #---------------------------------------- make a list if restarts.pbyd > 1
      newV = xData %*% W.list[[k]] # To be more accurate, we should update the index of W.list if restarts.pbyd > 1 
      lagV = newV
      for (i in 1:500) { # ---------- loop for relax_logistic
        iteration = i
        # update U:
        txv = t(xData)%*%newV
        svd.txv = La.svd(txv)
        newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U
        # update V:
        newV = xData%*%newU
        newV = prox_newton(v = newV, nu = nu)# here use the prox_newton function. 
        
        converge[iteration] = norm(lagV-newV) / ncol(lagV) # frobICA(S1=newV, S2=lagV) # BRISK note: normalize by number of columns
        lagV = newV
        if (converge[iteration] < eps) {break}
      } #--------------------------- end the loop for relax_logistic
      
      if (k != 1) {out.list[[k]] = out.list[[1]]}
      out.list[[k]]$loglik = sum(-lagV/(sqrt(3)/pi) - 2*log(1+exp(-lagV/(sqrt(3)/pi))))
      loglik[k] = out.list[[k]]$loglik
      out.list[[k]]$estS = newV
      out.list[[k]]$converge = converge
    } #----------------------------------------------------------- end the loop to make a list
    

    out.list = out.list[[which.max(loglik)]]
    out.list$distribution = "logistic"
    out.list$whitener = whitener
    
    if (converge_plot == F) {
      out.list = out.list[-3]
    } else {
      converge = out.list$converge
      plot = recordPlot()
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
      out.list$converge = plot
    }
    
  }
  
  ########################### 3
  if (method == "FastICA_logistic") {
    estX_logis = mlcaFP(xData = xData, n.comp = n.comp, whiten = "none", restarts.pbyd = restarts.pbyd, distribution='logistic', eps = eps)
    out.list$loglik = estX_logis$loglik
    out.list$estS = estX_logis$S
    out.list$distribution = estX_logis$distribution
    out.list$whitener = whitener
    if (converge_plot == T) {
      out.list$converge = "No converge plot for FastICA_logistic and FastICA_tanh"
    }
  }
  
  
  ########################### 4
  if (method == "FastICA_tanh") {
    estX_tanh = mlcaFP(xData = xData, n.comp = n.comp, whiten = "none", restarts.pbyd = restarts.pbyd, distribution='tanh', eps = eps)
    out.list$loglik = estX_tanh$loglik
    out.list$estS = estX_tanh$S
    out.list$distribution = estX_tanh$distribution
    out.list$whitener = whitener
    if (converge_plot == T) {
      out.list$converge = "No converge plot for FastICA_logistic and FastICA_tanh"
    }
  }
  
  return(out.list)
  
} # end of refinedICA function


