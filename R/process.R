#' Change the sign of S and M matrices to positive skewness.
#' @param S The S matrix with dimension P x Q.
#' @param M The M matrix with dimension Q x T.
#'
#' @return  A list of S and M matrices with positive skewness.
#' @export
signchange = function(S,M=NULL) {
  signskew = sign(apply(S,1,function(x) mean(x^3)))
  newS = signskew*S
  if(!is.null(M)) {
    newM = t(signskew*t(M))
  }else {newM=M}
  return(list(S=newS,M=newM))
}

#' For a given angle theta, returns a d x d Givens rotation matrix
#' @param theta The value of theta.
#' @param d The value of d.
#' @param which The value of which.
givens.rotation <- function(theta=0, d=2, which=c(1,2)){
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
  return(M)
}


#' Convert angle vector into orthodox matrix
#'
#' @param theta A vector of angles theta.
#' @return An orthodox matrix.
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
  return(W)
}




