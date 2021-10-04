#----------------------------------
# Benjamin Risk
# brisk@emory.edu
# Functions supporting LCA
# This is beta software. R package to come.
#------------------------------
mlcaFP <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB','tanh'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), max.comp = FALSE, reinit.max.comp=FALSE, df=0, irlba=FALSE,...) {

    #former option:
    #alg.typ = c('parallel','deflation'),
    #alg.typ = match.arg(alg.typ)
    alg.typ = 'parallel'

  
    distribution = match.arg(distribution)
    whiten=match.arg(whiten)

    if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
    ## whiten:
    
    if(irlba) require(irlba)
    if(max.comp) { #if statement evaluates to true for all max.comp!=0
      s.comp = n.comp
      n.comp = max.comp
    }
    if(max.comp=='TRUE') stop('max.comp should be an integer or FALSE')
    if(reinit.max.comp && max.comp==FALSE) stop('Can not reinitialize from max.comp solution if max.comp==FALSE')
    if(reinit.max.comp && alg.typ=='deflation') stop('reinit.max.comp not yet written for deflation algorithm')
    require(ProDenICA)
    #require(multidcov)
    if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
    if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
    if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')
    if(distribution=='logistic') Gfunc = logistic
    if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')
    if(distribution=='tanh'  && df>0) stop('df should be set to zero when using tanh')
    if(distribution=='tanh') Gfunc = gtanh
    
    if(!is.null(W.list) & class(W.list)!='list') stop('W.list must be a list')
    if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')

    orth.method= match.arg(orth.method)
    p = ncol(xData)
    nRow = nrow(xData)
    d = n.comp
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
    # warning('TO DO: Look at whitening methods and check for inconsistent options')
  if (is.null(W.list)) {
    if(restarts.pbyd) W.list = gen.inits(p=p,d=d,runs=restarts.pbyd,orth.method=orth.method)
    if(restarts.dbyd) {
      W.temp = gen.inits(p=d,d=d,runs=restarts.dbyd,orth.method=orth.method)
      #pad with zeros:
      zeros = matrix(0,p-d,d)
      W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros))
      W.list = c(W.list,W.temp)
    }
    
  }
  ## If restarts.pbyd and restarts.dbyd both equal zero:
  if (is.null(W.list)) W.list = gen.inits(p=p,d=d,runs=1,orth.method=orth.method)
  runs = length(W.list)
  out.list = NULL
  loglik.v = numeric(runs)
  for(k in 1:runs) {
    W0 = as.matrix(W.list[[k]])
    if(alg.typ == 'parallel') {
      out.list[[k]] = lca.par(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df, ...)
    out.list[[k]]$df = df
    }
    if(alg.typ == 'deflation') {
      out.list[[k]] = lca.def(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
    out.list[[k]]$df = df
    }
    if(max.comp) {
      flist0 = list()
      for (j in 1:d) flist0[[j]] <- Gfunc(out.list[[k]]$S[, j], ...)
      loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
      orderedLL = order(loglik.S,decreasing=TRUE)
      out.list[[k]]$S = out.list[[k]]$S[,orderedLL[1:s.comp]]
      out.list[[k]]$Ws = out.list[[k]]$Ws[,orderedLL[1:s.comp]]
      out.list[[k]]$loglik = sum(loglik.S[orderedLL[1:s.comp]])
      loglik.v[k] = out.list[[k]]$loglik
      } else {
      loglik.v[k] = out.list[[k]]$loglik
    }
  }
  for(i in 1:runs){
    out.list[[i]]$distribution=distribution
    out.list[[i]]$whitener = whitener
  }
  out = out.list[[which.max(loglik.v)]]
  if(reinit.max.comp) {
    out = lca.par(xData=xData,W0=out$Ws,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=s.comp,df=df, ...)
    out$df = df
    out.list[[k+1]] = out
  }
  if(out.all==TRUE) out.list else out
}
#---------------------------------

#-------------------------
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

#-------------------------------------

#-----------------------
# logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
#   #maximizes likelihood given s then calculates gradient w.r.t. w.hat
#   #df is not used
#   xData = as.vector(xData)
#   list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = (2*exp(-2*xData/scale) - 2*exp(-xData/scale)*(1+exp(-xData/scale))) / (scale^2*(1+exp(-xData/scale))^2))
# }

logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
  #maximizes likelihood given s then calculates gradient w.r.t. w.hat
  #df is not used
  xData = as.vector(xData)
  list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = - 2*exp(-xData/scale) / (scale^2*(1+exp(-xData/scale))^2))
}

# jb.stat written by Ze Jin:
jb.stat <- function(x, df=0) {
  n <- length(x)
  s <- sum(x^3)
  k <- sum(x^4)
  Gs <- x^3 * s / n^2 + (x^4 * k / n^2 + 9 / n - 6 * x^4 / n) / 4
  gs <- 6 * x^2 * s / n^2 + (8 * x^3 * (k / n - 3) / n) / 4
  gps <- 6 * (3 * x^4 + 2 * x * s) / n^2 + (24 * x^2 * (k / n - 3) / n + 32 * x^6 / n^2) / 4
  list(Gs = Gs, gs = gs, gps = gps)
}
 
gtanh <- function(s, df=0) {
  a=-1
  #warning('issue with signs -- not equivalent to fastICA, uses negative of G1 from ProDenICA')
  list(Gs = logb(cosh(a * s))/a, gs = tanh(a * s), gps = a * (1 - 
  tanh(a * s)^2))
}
 

#-----------------------
orthogonalize = function (W) {
  ##For arbitrary W, this is equivalent to (WW^T)^{-1/2} W
  temp <- svd(W)
  tcrossprod(temp$u,temp$v)
}

#-------------
orthonormalize <- function(xk,X,k) {
  #Gram-Schmidt orthogonalization
  #assumes columns of X have length equal to 1
  if(k!=1) {
    t <- numeric(length(xk))
    for (u in 1:(k-1)) {
      a <- sum(xk * X[,u])
      t <- t + a * X[,u]
    }
    xk <- xk - t
    }
  xk / sqrt(crossprod(xk))
}

genWfromWs <- function(Ws) {
  d = ncol(Ws)
  p = nrow(Ws)
  tempW = cbind(Ws,diag(p)[,(d+1):p])
  for(k in (d+1):p) {
    oldWk = tempW[,k]
    tempWk = tempW[,k]
    for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj,tempWj)
  }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  }
  tempW
}

temp.orthogonalize <- function(V,W) {
  #orthogonalizes the vector V to all columns in W
  #and returns cbind(W,orthV)
  oldWk = V
  tempWk = V
  tempW=cbind(W,V)
  k=ncol(W)+1
  for(j in 1:(k-1)) {
      tempWj = tempW[,j]
      tempWk = tempWk - tempWj * crossprod(tempWj,oldWk)/crossprod(tempWj)
    }
  tempW[,k] = tempWk/sqrt(crossprod(tempWk))
  tempW
}


#------------------------------------------------
# Parallel algorithm:
lca.par <- function(xData,W0,Gfunc,maxit,verbose,density,eps,n.comp,df=0,...) {
  W0 = as.matrix(W0)
  d = ncol(W0)
  if(n.comp!=d) stop('W0 needs to be p x d')
  p = ncol(xData)
  nRow = nrow(xData)
  s <- xData %*% W0
  flist <- as.list(1:d)
  ##  Densities of likelihood components:
  for (j in 1:d) flist[[j]] <- Gfunc(s[, j], df=df,...)
  
  
  flist0 <- flist
  crit0 <- mean(sapply(flist0, "[[", "Gs"))
  nit <- 0
  nw <- 10
  repeat {
    nit <- nit + 1
    gS <- sapply(flist0, "[[", "gs")
    gpS <- sapply(flist0, "[[", "gps")
    #t1 <- t(xData) %*% gS/nRow
    t1 <- crossprod(xData,gS)/nRow
    t2 <- apply(gpS, 2, mean)
    if(d>1) W1 <- t1 - W0%*%diag(t2) else W1 <- t1 - W0*t2
    W1 <- orthogonalize(W1)
    if(d>1) nw <- frobICA(t(W0), t(W1))^2 else nw <- mean((W0-W1)^2) #Uses a measure that works for non-square matrices -- MSE. The measure is defined for M so here we use transpose of W.
    W0 <- W1
    s <- xData %*% W0
    for (j in 1:d) flist0[[j]] <- Gfunc(s[, j], df=df, ...)
    crit0 <- mean(sapply(flist0, "[[", "Gs"))
    if (verbose) cat("Iter", nit, "G", crit0, "Delta", nw, "\n")
    if ((nit >= maxit)) {
      warning('Max iter reached')
      break
    }
    if (nw < eps) break
  }
  out = list(Ws = W0, loglik = d*nRow*crit0, S = s)
  if(density) out$density = lapply(flist0, "[[", "density")
  out
}

#--------------------------------------
myMixmat <-  function (p = 2) {
  a <- matrix(rnorm(p * p), p, p)
  sa <- svd(a)
  d <- sort(runif(p,min=1,max=10))
  mat <- sa$u %*% (sa$v * d)
  attr(mat, "condition") <- d[p]/d[1]
  mat
}

#------------------------
#---------------------------------------
standardizeM <- function(M) {
  #M is d x p
  diag(diag(M%*%t(M))^(-1/2))%*%M
}

#--------------------
tiltedgaussian = function (xData, df = 8, B = 100, ...) {
  #This function is based on ProDenICA::GPois by Trevor Hastie
  #NOTE: Assumes data are zero mean.
  require(gam)
  n <- length(xData)
  sd.x = sd(xData)
  rx <- c(min(xData)-0.1*sd.x, max(xData)+0.1*sd.x) 
  xg <- seq(from = rx[1], to = rx[2], length = B)
  gaps <- diff(rx)/(B - 1)
  xcuts <- c(rx[1] - gaps/2, xg[-B] + gaps/2, rx[2] + gaps/2)
  #NOTE: I use the response variable that corresponds to the LCA paper.
  #This differs from the GPois algorithm in ProDenICA
  ys <- as.vector(table(cut(xData, xcuts)))/(gaps*n)
  pois.fit <- suppressWarnings(gam(ys ~ s(xg, df)+offset(dnorm(xg,log=TRUE)), family = poisson, ...))
  Gs <- predict(pois.fit) #log tilt function predicted at grid locations (note: predict on gam object can not be used to obtain derivatives)
  # the gam object with the predict function can not be used directly to obtain the derivatives
  # of the smoothing spline.
  # Here, we refit another iteration of the IRWLS algorithm used in gam:
  # Note: working residuals = (y - mu0)/mu0
  # weights = mu0
  # fitted(pois.fit) = mu0
  # predict(pois.fit) = eta0 = log(mu0)
  sGs = Gs #+ log(sum(dnorm(xg))/sum(fitted(pois.fit)))
  z0 <- sGs + residuals(pois.fit, type='working')
  pois.refit <- smooth.spline(x=xg, y=z0, w=fitted(pois.fit),df=df) #obtain the log tilt function in an object that can be used to obtain derivatives
  Gs <- predict(pois.refit, xData, deriv = 0)$y
  gs <- predict(pois.refit, xData, deriv = 1)$y
  gps <- predict(pois.refit, xData, deriv = 2)$y
  fGs <- function(x) predict(pois.refit,x,deriv=0)$y
  fgs <- function(x) predict(pois.refit,x,deriv=1)$y
  fgps <- function(x) predict(pois.refit,x,deriv=2)$y
  list(Gs = Gs, gs = gs, gps = gps, fGs = fGs, fgs=fgs, fgps=fgps)
}
#---------------------------------------------

#----------------------------
# estimate mixing matrix from estimates of components:
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}


#-----------------------------------
# order by likelihood
# option for positive skewness
order.likelihood <- function(S,positive.skew=TRUE,distribution=c('logistic','tiltedgaussian','logcosh'),out.loglik=FALSE,...) {
  distribution = match.arg(distribution)
  nObs = nrow(S)
  d = ncol(S)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(positive.skew) {
    skewness <- function(x, n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    skew = apply(S, 2, skewness)
    sign = -1 * (skew < 0) + 1 * (skew > 0)
    S = S %*% diag(sign)
  }
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
  orderedLL = order(loglik.S,decreasing=TRUE)
  S = S[,orderedLL]
  if(out.loglik) return(list(S=S,loglik=sort(loglik.S,decreasing=TRUE))) else S
}

marginal.likelihoods <- function(S,distribution=c('logistic','tiltedgaussian','logcosh','GPois'),...)
{
  distribution = match.arg(distribution)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='logcosh') Gfunc = ProDenICA::G1
  if(distribution=='GPois') Gfunc = ProDenICA::GPois
  d = ncol(S)
  flist0 = list()
  for (j in 1:d) flist0[[j]] <- Gfunc(S[, j], ...)
  apply(sapply(flist0, "[[", "Gs"),2,sum)
}



#-------------------------
#Match based on L2 distances
matchICA=function(S,template,M=NULL) {
  require(clue)
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

angleMatchICA=function(Mx,My,Sx=NULL,Sy=NULL) {
  # match the colums of Mx and Mx, using the 
  # n x p parameterization of the JIN decomposition
  require(clue)
  rx=ncol(Mx)
  ry=ncol(My)
  n.comp = max(rx,ry)
  
  n.obs=nrow(Mx)
  if(n.comp>n.obs) warning('Input should be n x r')
  if(n.obs!=nrow(My)) warning('Mx and My need to have the same number of rows')
  Mx = t(Mx)
  Mx = Mx / sqrt(apply(Mx^2,1,sum))
  
  
  My = t(My)
  My = My / sqrt(apply(My^2,1,sum))
  
  if(rx<ry) {
    Mx = rbind(Mx,matrix(0,ry-rx,n.obs))
  }
  if(ry<rx) {
    My = rbind(My,matrix(0,rx-ry,n.obs))
  }
  angle.mat1=acos(Mx%*%t(My))
  angle.mat2=acos(-1*Mx%*%t(My))
  angle.mat=angle.mat1*(angle.mat1<=angle.mat2)+angle.mat2*(angle.mat2<angle.mat1)
  map=as.vector(solve_LSAP(angle.mat))
  angle.mat1.perm = angle.mat1[,map]
  angle.mat2.perm = angle.mat2[,map]
  angle1=diag(angle.mat1.perm)
  angle2=diag(angle.mat2.perm)
  matchedangles = apply(cbind(angle1,angle2),1,min)
  allangles = angle.mat1.perm*(angle.mat1.perm<=angle.mat2.perm)+angle.mat2.perm*(angle.mat2.perm<angle.mat1.perm)
  sign.change=-1*(angle2<angle1)+1*(angle1<=angle2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)
  
  My.perm=t(perm)%*%My

  # reorder components by their matched angles 
  smatchedangles = sort(matchedangles)
  omangles = order(matchedangles)
  sallangles = allangles[omangles[1:rx],omangles[1:ry]]
  sMx = Mx[omangles[1:rx],]
  sMy.perm = My.perm[omangles[1:ry],]

  if(!is.null(Sy)) {
  Sy.perm=t(perm)%*%Sy
  sSy.perm = Sy.perm[omangles[1:ry],]
  sSx = Sx[omangles[1:rx],]
  return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles,Sx = sSx, Sy = sSy.perm))
  }
  else {
    return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles))
  }
}






#-------------------------------------
# Match mixing matrices:
# This function does not require M to be square:
frobICA<-function(M1=NULL,M2=NULL,S1=NULL,S2=NULL,standardize=FALSE) {
  #MODEL: X = S M + E, so M is d x p
  #standardize: if standardize==TRUE, then standardizes rows of M1 and M2
  #to have unit norm; if using S1 and S2, standardizes columns to have unit variance.
  #standardize=TRUE makes the measure scale invariant.

  require(clue)
  tfun = function(x) all(x==0)
  if(is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) stop("need to supply either M1 and M2 or S1 and S2")
  if(!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if(!is.null(M1) && nrow(M1) > ncol(M1)) stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")

  if(is.null(M1)) {
    nS = nrow(S1)
    if(nS!=nrow(S2)) stop('S1 and S2 must have the same number of rows')
    if(sum(apply(S1,2,tfun)) + sum(apply(S2,2,tfun))) stop('frobICA not defined when S1 or S2 has a column of all zeros')
    if(standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if(p < q) {
      S1 = cbind(S1,matrix(0,nS,(q-p)))
    }
    if(q < p) {
      S2 = cbind(S2,matrix(0,nS,(p-q)))
    }
    Stemp = matchICA(S=S1,template=S2)
    n.comp = max(q,p)
    indices = c(1:n.comp)[!(apply(Stemp,2,tfun) | apply(S2,2,tfun))]
    return(sqrt(sum((Stemp[,indices] - S2[,indices])^2))/sqrt(nS*min(p,q)))
  }

  else {
    if(sum(apply(M1,1,tfun)) + sum(apply(M2,1,tfun))) stop('frobICA not defined when M1 or M2 has a row of all zeros')
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
    p = ncol(M1)
    if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp=max(d,q)
    if(n.comp > p) warning("M should be d x p")
    if(d<q) {
      M1 = rbind(M1,matrix(0,(q-d),p))
    }
    if(q<d) {
      M2 = rbind(M2,matrix(0,(d-q),p))
    }
    l2.mat1=l2.mat2=matrix(NA,nrow=n.comp,ncol=n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        #since signs are arbitrary, take min of plus and minus:
        l2.mat1[i,j]=sum((M2[i,]-M1[j,])^2)
        l2.mat2[i,j]=sum((M2[i,]+M1[j,])^2)
      }
    }
    l2.mat1=sqrt(l2.mat1)
    l2.mat2=sqrt(l2.mat2)
    #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
    l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
    map=as.vector(solve_LSAP(l2.mat))
    #retain relevant l2 distances:
    l2.1=diag(l2.mat1[,map])
    l2.2=diag(l2.mat2[,map])
    #sign.change is for re-ordered matrix 2
    sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
    perm=diag(n.comp)[,map]%*%diag(sign.change)
    M.perm=t(perm)%*%M1
    indices = c(1:n.comp)[!(apply(M.perm,1,tfun) | apply(M2,1,tfun))]
    return(sqrt(sum((M.perm[indices,]-M2[indices,])^2))/sqrt(p*min(d,q)))
  }
}


#----------------
#----------------------------------------
#Function to make most extreme values for the skewed tail positive, i.e., force all distributions to be right skewed, and order ICs by skewness.
rightskew=function(S,M=NULL,order.skew=TRUE) {
  #S: n x d matrix
  #A: d x d` corresponding to X = S A
  #If order = TRUE, then ICs are organized from HIGHEST to LOWEST skewness where skewness is forced to be positive for all ICs.
  nObs <- nrow(S)
  if(ncol(S)>nObs) stop('S must be n x d')
  skewness<-function(x,n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  skew=apply(S,2,skewness)
  sign=-1*(skew<0)+1*(skew>0)
  S.new=S%*%diag(sign)
  if(!is.null(M)) M.new=t(diag(sign))%*%M
  if(order.skew==TRUE) {
    skew.new=apply(S.new,2,skewness)
    perm=order(skew.new,decreasing=TRUE)
    S.new=S.new[,perm]
    if(!is.null(M)) M.new=M.new[perm,]
  }
  if(is.null(M)) {
    S.new
  } else
    return(list(S=S.new,M=M.new))
}
#---------------------------------------------

mProDenICA <- function(X, n.comp = ncol(X), restarts=0, tol=1e-07,maxit=100,G = c('GPois','G0','G1'),verbose=FALSE,whiten=FALSE,...) {
  ##NOTE: the restarts in ProDenICA evaluate the likelihood at a sample of orthogonal matrices, identifies the random matrix associated with highest likelihood, and then estimates ICs for this single initialization. Here, I initiate from the entire set of random matrices.
  ##NOTE: Restarts defined differently here than in ProDenICA. ProDenICA is initiatialized from restarts+1 initial values.
  ##NOTE: G defined differently from ProDenICA's Gfunc; here it is a string
  require(ProDenICA)
  G = match.arg(G)
  if(G=='G0') Gfunc=G0
  if(G=='G1') Gfunc=G1
  if(G=='GPois') Gfunc=GPois
  est.list = list()
  runs = restarts+1
  obj.v = numeric(runs)
  theta.list = lapply(rep(choose(n.comp, 2), runs), runif, min = 0, max = 2 * pi)
  W.list = lapply(theta.list, theta2W)
  if(whiten) {
    a<- whitener(X=X,n.comp=n.comp)
    zData <- a$Z
  } else {
    zData <- X[,1:n.comp]
  }

  for(i in 1:runs) {
    est.list[[i]] = ProDenICA(x=zData, k=n.comp, W0=W.list[[i]], whiten=FALSE, maxit = maxit, thresh = tol, trace=verbose, restarts=0, Gfunc=Gfunc,...)
    if (G=='G1') obj.v[i] = calc.negent.hyvarinen(s=est.list[[i]]$s) else obj.v[i] = est.list[[i]]$negentropy
  }
  out = est.list[[which.max(obj.v)]]
  if(G=='G1') out$negentropy=obj.v[which.max(obj.v)]
  if(whiten) out$whitener=a$whitener
  out
}

theta2W = function(theta)
{
  #<author hidden>
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

# Whitening Function:
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    require(irlba)
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}


# Edits 10 March 2017
# Returns square root of the precision matrix:
covwhitener <- function(X,n.comp=ncol(X),center.row=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Creates model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  covmat = cov(x.center)
  evdcov = eigen(covmat,symmetric = TRUE)
  whitener = evdcov$vectors%*%diag(1/sqrt(evdcov$values))%*%t(evdcov$vectors)
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=whitener,Z=x.center%*%whitener,mean=apply(X,2,mean)))
}



rtwonorm <- function(n, mean=c(0,5), sd=c(2/3,1), prob=0.5) {
  k <- rbinom(n,1,prob=prob)
  k*rnorm(n,mean[1],sd[1])+(1-k)*rnorm(n,mean[2],sd[2])
}

rmixnorm <- function(n, pars = list(mean=c(0,5), sd = c(2/3,1), prob=c(0.25,0.75))) {
  probs = pars[['prob']]
  means = pars[['mean']]
  sigmas = pars[['sd']]
  if(sum(probs)!=1) stop('Probabilities must sum to one')
  z = rmultinom(n=n, size=1, prob = probs)
  k = length(probs)
  #use rnorm recycling:
  x = rnorm(k*n,means,sigmas)
  dim(x) = c(k,n)
  apply(x*z,2,sum)
}

SimTwoNorms <- function(n.samples, distribution=c('mix-sub','mix-super'),snr,noisyICA=FALSE) {
  distribution = match.arg(distribution)
  if(distribution=='mix-sub') {
    mean=c(-1.7,1.7); sd=c(1,1); prob=0.75
  }
  if(distribution=='mix-super') {
    mean=c(0,5); sd=c(2/3,1); prob=0.95
  }
  sim.S <- rtwonorm(n=2*n.samples, mean=mean, sd=sd, prob=prob)
  dim(sim.S) <- c(n.samples,2)
  sim.M = myMixmat(5)
  sim.W = solve(sim.M)
  if(noisyICA) {
    sim.N <- matrix(rnorm(n=5*n.samples,mean=0,sd=1),nrow=n.samples,ncol=5)
  } else {
    sim.N <- matrix(rnorm(n=3*n.samples,mean=0,sd=1),nrow=n.samples,ncol=3)
  }

  sim.Ms = sim.M[1:2,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn <- sim.N
  } else {
    sim.Mn <- sim.M[3:5,]
    sim.Xn <- sim.N%*%sim.Mn
  }
  #alpha = 1/sqrt(mean(sim.Xs^2))
  alpha = 1/sd(as.vector(sim.Xs))
  sim.Xs = sim.Xs*alpha
  mean.S = apply(sim.S,2,mean)
  temp.S = scale(sim.S,center=TRUE,scale=FALSE)
  scalingMat = diag(apply(temp.S,2,sd))
  scaled.sim.S = temp.S%*%solve(scalingMat)
  scaled.sim.Ms = sqrt(snr)*alpha*scalingMat%*%sim.Ms
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn))
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S),scaled.sim.Ms = scaled.sim.Ms,scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}


SimLCA <- function(n.samples, distribution = c('logistic','t','gumbel'), nu.vector = c(sqrt(3)/pi,sqrt(3)/pi), snr, k = NULL, noisyICA=FALSE) {
  ##k = number of noise components
  require(MASS)
  if(is.null(k)) {
    if (noisyICA) k = 5 else k = 3
  }
  d = length(nu.vector)
  p = ifelse(noisyICA==TRUE,k,k+d)
  distribution=match.arg(distribution)
  sim.S = NULL
  if(distribution == 'logistic') for (i in 1:d) sim.S = cbind(sim.S,rlogis(n=n.samples,scale=nu.vector[i]))
  if(distribution == 't') for (i in 1:d) sim.S = cbind(sim.S,rt(n=n.samples,df=nu.vector[i]))
  if(distribution == 'gumbel') {
    require(evd)
    for (i in 1:d) sim.S = cbind(sim.S,rgumbel(n=n.samples,scale=nu.vector[i]))
  }
  sim.N = matrix(rnorm(n.samples*k,sd=1),nrow=n.samples)
  sim.M = myMixmat(p)
  sim.Ms = sim.M[1:d,]
  sim.Xs = sim.S%*%sim.Ms
  if(noisyICA) {
    sim.Mn = NULL
    sim.Xn = sim.N
  } else {
    sim.Mn = sim.M[(d+1):p,]
    sim.Xn = sim.N%*%sim.Mn
  }
  #alpha = 1/sqrt(mean(sim.Xs^2))
  alpha = 1/sd(as.vector(sim.Xs))
  sim.Xs = sim.Xs*alpha
  mean.S = apply(sim.S,2,mean)
  temp.S = scale(sim.S,center=TRUE,scale=FALSE)
  scalingMat = diag(apply(temp.S,2,sd))
  scaled.sim.S = temp.S%*%solve(scalingMat)
  scaled.sim.Ms = sqrt(snr)*alpha*scalingMat%*%sim.Ms
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2))
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) #eigenvalues sum to p
  sim.Xs = sqrt(snr)*sim.Xs #equivalent to scaled.sim.S%*%(alpha*sqrt(snr)*scalingMat%*%sim.Ms)+alpha*sqrt(snr)*matrix(mean.S,nrow=n.samples,ncol=2,byrow=TRUE)%*%sim.Ms
  #since we only recover scaled.sim.S, "true Ms" for, e.g., IFA, is defined as in scaled.sim.Ms
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)

  return(list(sim.S = sim.S, sim.N = sim.N, sim.Ms = sim.Ms, scaled.sim.Ms = scaled.sim.Ms, sim.Mn = sim.Mn, sim.X=sim.X, scaled.sim.S = scale(sim.S), scaled.sim.X = scale(sim.X), whitened.sim.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
}



################################################
SimFMRI123 = function(snr = 1, noisyICA=FALSE, nTR=50, nImages=1, phi=0.5, dim.data=c(33,33), var.inactive=0.0001) {
  ##ASSUME 1,000 samples
  require(neuRosim)
  require(steadyICA)
  m = nImages
  #Latent components are fixed for each simulation:
  x1 = rep(3,5)
  y1 = c(3:7)
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))
  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)
  
  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))
  
  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)
  
  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))
  
  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3))
  
  if(m>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }
  
  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S[sim.S==0] = baseline
  
  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:m) {
    t.sim.GRF <- spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    dim(t.sim.GRF) <- c(prod(dim.data),nscan)
    sim.GRF = rbind(sim.GRF,t.sim.GRF)
  }
  
  ##Mixmat:
  #create timecourses for latent components:
  totaltime <- nTR
  nOnsets = 5+1
  onsets <- seq(from=1, to=totaltime, length=nOnsets)
  dur <- totaltime/10
  #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  #NOTE: Time courses can not be identical.
  row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')
  
  sim.Ms = matrix(c(row1,row2,row3),nrow=3,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms
  
  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2)) 
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2)) 
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) 
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  
  if(noisyICA) { 
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}
#--------------------------------------------

# function to create network matrices from vectorized lower diagonals:

vec2net = function(invector) {
  #invector: p x 1, where p is the number of nodes
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  diag(outNet) = 1
  outNet
}


