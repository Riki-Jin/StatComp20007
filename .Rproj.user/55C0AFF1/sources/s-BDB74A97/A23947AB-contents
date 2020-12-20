#' @title A simplified Mardia multivariate skewness and kurtosis calculator
#' @description A simplified Mardia multivariate skewness and kurtosis calculator
#' @param x A data matrix
#' @param na.rm An indication of the missing data, the default value is True
#' @return Mardia multivariate skewness and kurtosis
#' @importFrom stats na.omit cov
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib StatComp20007
#' @examples 
#' \dontrun{
#' x <- mvrnorm(1000,c(1,2,3),matrix(c(4,0.7,0.3,0.7,2,1.4,0.3,1.4,5),nrow=3,ncol=3))
#' mardias(x)
#' }
#' @export

mardias <- function(x, na.rm = TRUE){    
  if (na.rm)
    x <- na.omit(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  x <- scale(x, scale = FALSE)
  S <- cov(x)*(n-1)/n
  S.inv <- MASS::ginv(S)
  D <- x %*% S.inv %*% t(x)
  b1p <- sum(D^3)/n^2
  b2p <- sum(diag(D^2))/n
  return(c(b1p,b2p))
}

#' @title Moment calculator of multivariate normal distribution
#' @description Moment calculator of multivariate normal distribution
#' @param p Dimension
#' @param sig p*p multivariate normal covariance matrix 
#' @return Moments flattened into an vector
#' @import symmoments
#' @examples 
#' \dontrun{
#' p <- 2
#' sig <- matrix(c(2,0.7,0.7,1),2,2)
#' A <- mstore(p,sig)
#' matrix(A,13,13)
#' }
#' @export 

mstore <- function(p,sig) {
  siglower <- sig[lower.tri(sig,diag=TRUE)]
  decito13 <- function(n,p){
    result <- rep(0,p)
    temp <- n
    for (i in 1:(p-1)){
      result[i] <- temp %/% (13^(p-i))
      temp <- temp-(13^(p-1)*result[i])
    }
    result[p] <- n %% 13
    return(result)
  }
  todeci <- function(a){
    p <- length(a)
    f <- rep(0,p)
    for (i in 1:p) 
      f[i] <- 13^(p-i)
    return(sum(a*f))
  }
  Astore <- rep(0,13^p)
  for (i in 0:(13^p-1)) {
    comb <- decito13(i,p)
    combsum <- sum(comb)
    if ((combsum<=12) & (combsum %% 2 == 0))
      Astore[i+1] <- evaluate(callmultmoments(comb),siglower)
  }
  Astore
}

#' @title Skewness calculator given Gaussian copula and Fleishman parameters
#' @description Skewness calculator given Gaussian copula and Fleishman parameters
#' @param p Dimension
#' @param sig p*p multivariate normal covariance matrix
#' @param b Fleishman parameter
#' @param c Fleishman parameter
#' @param d Fleishman parameter
#' @return Exact mardia multivariate skewness
#' @import symmoments
#' @import MASS
#' @examples 
#' \dontrun{
#' p <- 2
#' sig <- matrix(c(1,0.7,0.7,1),2,2)
#' b <- 1
#' c <- 2
#' d <- 3
#' gauss_skew(p,sig,b,c,d)
#' }
#' @export 

gauss_skew <- function(p,sig,b,c,d) {
  decito13 <- function(n,p){
    result <- rep(0,p)
    temp <- n
    for (i in 1:(p-1)){
      result[i] <- temp %/% (13^(p-i))
      temp <- temp-(13^(p-1)*result[i])
    }
    result[p] <- n %% 13
    return(result)
  }
  
  todeci <- function(a){
    p <- length(a)
    f <- rep(0,p)
    for (i in 1:p) 
      f[i] <- 13^(p-i)
    return(sum(a*f))
  }
  
  gauss_v <- function(p,sig,b,c,d) {
    S <- matrix(0,p,p)
    a <- -c    # SCALE
    abcd <- c(a,b,c,d)
    siglower <- sig[lower.tri(sig,diag=TRUE)]
    for (i in 1:p)
      for (j in 1:p)
        for (ii in 1:4)
          for (jj in 1:4)
            if ((ii+jj) %% 2 == 0) {
              xx <- rep(0,p)
              xx[i] <- xx[i]+(ii-1)
              xx[j] <- xx[j]+(jj-1)
              S[i,j] <- S[i,j] + abcd[ii]*abcd[jj]*
                evaluate(callmultmoments(xx),siglower)
            }
    V <- solve(S)
    V
  }
  
  Astore <- mstore(p,sig)
  A <- array(0,c(p,p,p))
  a <- -c
  abcd <- c(a,b,c,d)
  V <- gauss_v(p,sig,b,c,d)
  for (i in 1:p)
    for (j in 1:p)
      for (k in 1:p)
        if (identical( c(i,j,k), sort(c(i,j,k)) )) {
          for (ii in 1:4)
            for (jj in 1:4)
              for (kk in 1:4)
                if ((ii+jj+kk) %% 2 == 1) {
                  xx <- rep(0,p)
                  xx[i] <- xx[i]+(ii-1)
                  xx[j] <- xx[j]+(jj-1)
                  xx[k] <- xx[k]+(kk-1)
                  A[i,j,k] <- A[i,j,k] + abcd[ii]*abcd[jj]*abcd[kk]*Astore[todeci(xx)+1]
                }
        } else {
          st <- sort(c(i,j,k))
          A[i,j,k] = A[st[1],st[2],st[3]]
        }
  skew <- 0
  for (i in 1:p)
    for (j in 1:p)
      for (i1 in 1:p)
        for (j1 in 1:p)
          for (i2 in 1:p)
            for (j2 in 1:p)
              skew <- skew + V[i,j]*V[i1,j1]*V[i2,j2]* 
    A[i,i1,i2]*A[j,j1,j2]
  return(skew)
}

#' @title Kurtosis calculator given Gaussian copula and Fleishman parameters
#' @description Kurtosis calculator given Gaussian copula and Fleishman parameters
#' @param p Dimension
#' @param sig p*p multivariate normal covariance matrix
#' @param b Fleishman parameter
#' @param c Fleishman parameter
#' @param d Fleishman parameter
#' @return Exact mardia multivariate kurtosis
#' @import symmoments
#' @import MASS
#' @examples 
#' \dontrun{
#' p <- 2
#' sig <- matrix(c(1,0.7,0.7,1),2,2)
#' b <- 1
#' c <- 2
#' d <- 3
#' gauss_kurt(p,sig,b,c,d)
#' }
#' @export 

gauss_kurt <- function(p,sig,b,c,d) {
  decito13 <- function(n,p){
    result <- rep(0,p)
    temp <- n
    for (i in 1:(p-1)){
      result[i] <- temp %/% (13^(p-i))
      temp <- temp-(13^(p-1)*result[i])
    }
    result[p] <- n %% 13
    return(result)
  }
  
  todeci <- function(a){
    p <- length(a)
    f <- rep(0,p)
    for (i in 1:p) 
      f[i] <- 13^(p-i)
    return(sum(a*f))
  }
  
  gauss_v <- function(p,sig,b,c,d) {
    S <- matrix(0,p,p)
    a <- -c    # SCALE
    abcd <- c(a,b,c,d)
    siglower <- sig[lower.tri(sig,diag=TRUE)]
    for (i in 1:p)
      for (j in 1:p)
        for (ii in 1:4)
          for (jj in 1:4)
            if ((ii+jj) %% 2 == 0) {
              xx <- rep(0,p)
              xx[i] <- xx[i]+(ii-1)
              xx[j] <- xx[j]+(jj-1)
              S[i,j] <- S[i,j] + abcd[ii]*abcd[jj]*
                evaluate(callmultmoments(xx),siglower)
            }
    V <- solve(S)
    V
  }
  
  Astore <- mstore(p,sig)
  A <- array(0,c(p,p,p,p))
  a <- -c
  abcd <- c(a,b,c,d)
  V <- gauss_v(p,sig,b,c,d)
  for (i in 1:p)
    for (j in 1:p)
      for (k in 1:p)
        for (l in 1:p)
          if (identical( c(i,j,k,l), sort(c(i,j,k,l)) )) {
            for (ii in 1:4)
              for (jj in 1:4)
                for (kk in 1:4)
                  for (ll in 1:4)
                    if ((ii+jj+kk+ll) %% 2 == 0) {
                      xx <- rep(0,p)
                      xx[i] <- xx[i]+(ii-1)
                      xx[j] <- xx[j]+(jj-1)
                      xx[k] <- xx[k]+(kk-1)
                      xx[l] <- xx[l]+(ll-1)
                      A[i,j,k,l] <- A[i,j,k,l] + abcd[ii]*abcd[jj]*abcd[kk]*abcd[ll]*Astore[todeci(xx)+1]
                    }
          } else {
            st <- sort(c(i,j,k,l))
            A[i,j,k,l] = A[st[1],st[2],st[3],st[4]]
          }
  kurt <- 0
  for (i in 1:p)
    for (j in 1:p)
      for (i1 in 1:p)
        for (j1 in 1:p)
          kurt <- kurt + V[i,j]*V[i1,j1]*A[i,i1,j,j1]
  return(kurt)
}

#' @title Non-normal random numbers generator with desired Mardia skewness and kurtosis.
#' @description Non-normal random numbers generator with desired Mardia skewness and kurtosis.
#' @param p Dimension
#' @param n Sample size
#' @param sig p*p multivariate normal covariance matrix
#' @param Sigma p*p desired result covariance matrix
#' @param ms Desired Mardia skewness
#' @param mk Desired Mardia kurtosis
#' @param initial A vector with 3 numbers for initial polynominal coefficients' (b,c,d). The default setting is (0.9,0.4,0).
#' @return A data matrix (non-normal random numbers generated)
#' @import symmoments
#' @import MASS
#' @import expm
#' @import stats
#' @examples 
#' \dontrun{
#' set.seed(1)
#' p <- 2
#' n <- 10000
#' ms <- 3
#' mk <- 20
#' sig <- matrix(c(1,0.7,0.7,1),2,2)
#' Sigma <- matrix(c(1,0.5,0.5,1),2,2)
#' }
#' @export

gaussnonr <- function(p,n,ms,mk,sig,Sigma,initial=NULL) {
  
  if (is.null(initial))
    initial <- c(0.9, 0.4, 0)
  siglower <- sig[lower.tri(sig,diag=TRUE)]
  
  decito13 <- function(n,p){
    result <- rep(0,p)
    temp <- n
    for (i in 1:(p-1)){
      result[i] <- temp %/% (13^(p-i))
      temp <- temp-(13^(p-1)*result[i])
    }
    result[p] <- n %% 13
    return(result)
  }
  
  todeci <- function(a){
    p <- length(a)
    f <- rep(0,p)
    for (i in 1:p) 
      f[i] <- 13^(p-i)
    return(sum(a*f))
  }
  
  Astore <- rep(0,13^p)
  for (i in 0:(13^p-1)) {
    comb <- decito13(i,p)
    combsum <- sum(comb)
    if ((combsum<=12) & (combsum %% 2 == 0))
      Astore[i+1] <- evaluate(callmultmoments(comb),siglower)
  }
  
  gauss_v <- function(p,sig,b,c,d) {
    S <- matrix(0,p,p)
    a <- -c    # SCALE
    abcd <- c(a,b,c,d)
    siglower <- sig[lower.tri(sig,diag=TRUE)]
    for (i in 1:p)
      for (j in 1:p)
        for (ii in 1:4)
          for (jj in 1:4)
            if ((ii+jj) %% 2 == 0) {
              xx <- rep(0,p)
              xx[i] <- xx[i]+(ii-1)
              xx[j] <- xx[j]+(jj-1)
              S[i,j] <- S[i,j] + abcd[ii]*abcd[jj]*
                evaluate(callmultmoments(xx),siglower)
            }
    V <- solve(S)
    V
  }
  
  gauss_skew <- function(Astore,p,sig,b,c,d) {
    A <- array(0,c(p,p,p))
    a <- -c
    abcd <- c(a,b,c,d)
    V <- gauss_v(p,sig,b,c,d)
    for (i in 1:p)
      for (j in 1:p)
        for (k in 1:p)
          if (identical( c(i,j,k), sort(c(i,j,k)) )) {
            for (ii in 1:4)
              for (jj in 1:4)
                for (kk in 1:4)
                  if ((ii+jj+kk) %% 2 == 1) {
                    xx <- rep(0,p)
                    xx[i] <- xx[i]+(ii-1)
                    xx[j] <- xx[j]+(jj-1)
                    xx[k] <- xx[k]+(kk-1)
                    A[i,j,k] <- A[i,j,k] + abcd[ii]*abcd[jj]*abcd[kk]*Astore[todeci(xx)+1]
                  }
          } else {
            st <- sort(c(i,j,k))
            A[i,j,k] = A[st[1],st[2],st[3]]
          }
    skew <- 0
    for (i in 1:p)
      for (j in 1:p)
        for (i1 in 1:p)
          for (j1 in 1:p)
            for (i2 in 1:p)
              for (j2 in 1:p)
                skew <- skew + V[i,j]*V[i1,j1]*V[i2,j2]* 
      A[i,i1,i2]*A[j,j1,j2]
    return(skew)
  }
  
  gauss_kurt <- function(Astore,p,sig,b,c,d) {
    A <- array(0,c(p,p,p,p))
    a <- -c
    abcd <- c(a,b,c,d)
    V <- gauss_v(p,sig,b,c,d)
    for (i in 1:p)
      for (j in 1:p)
        for (k in 1:p)
          for (l in 1:p)
            if (identical( c(i,j,k,l), sort(c(i,j,k,l)) )) {
              for (ii in 1:4)
                for (jj in 1:4)
                  for (kk in 1:4)
                    for (ll in 1:4)
                      if ((ii+jj+kk+ll) %% 2 == 0) {
                        xx <- rep(0,p)
                        xx[i] <- xx[i]+(ii-1)
                        xx[j] <- xx[j]+(jj-1)
                        xx[k] <- xx[k]+(kk-1)
                        xx[l] <- xx[l]+(ll-1)
                        A[i,j,k,l] <- A[i,j,k,l] + abcd[ii]*abcd[jj]*abcd[kk]*abcd[ll]*Astore[todeci(xx)+1]
                      }
            } else {
              st <- sort(c(i,j,k,l))
              A[i,j,k,l] = A[st[1],st[2],st[3],st[4]]
            }
    kurt <- 0
    for (i in 1:p)
      for (j in 1:p)
        for (i1 in 1:p)
          for (j1 in 1:p)
            kurt <- kurt + V[i,j]*V[i1,j1]*A[i,i1,j,j1]
    return(kurt)
  }
  
  optimtarget_gauss <- function(x, a){
    ##Revised Fleishman method
    b = x[1]
    c = x[2]
    d = x[3]
    g1 = a[[1]][1]
    g2 = a[[1]][2]
    p = a[[1]][3]
    sig = a[[2]]
    Astore <- a[[3]]
    z = 
      (1 - (b^2+2*c^2+6*b*d+15*d^2))^2 +
      (g1 - gauss_skew(Astore,p,sig,b,c,d))^2 +
      (g2 - gauss_kurt(Astore,p,sig,b,c,d))^2
  }
  
  findcoe_gauss = function(canshu,initial){
    ##
    #Uses the built in minimization function to solve for b, c, and d#
    # if the third and fourth moment of ksi are given.
    ##
    output = optim(initial,optimtarget_gauss,a = canshu,method = "BFGS",
                   control = list(ndeps = rep(1e-10, 3),reltol = 1e-10,maxit = 1e8))$par
    
    return(output)
  }
  
  a = list(c(ms,mk,p),sig,Astore)
  bcd = findcoe_gauss(a,initial)
  Z = mvrnorm(n,rep(0,p),sig)
  b = bcd[1]
  c = bcd[2]
  d = bcd[3]
  xi = -c+b*Z+c*Z^2+d*Z^3
  V = gauss_v(p,sig,b,c,d)
  V_h = sqrtm(V)
  eta = matrix(0, nrow = n, ncol = p)
  for (j in 1:n){
    for(m in 1:p){
      for(i in 1:p){
        eta[j,m] = eta[j,m] + V_h[i,m] * xi[j,i]}}}
  #  eta = V_h %*% t(xi)
  x = matrix(0, nrow = n, ncol = p)
  r = chol(Sigma)
  
  for (j in 1:n){
    for(m in 1:p){
      for(i in 1:p){
        x[j,m] = x[j,m] + r[i,m] * eta[j,i]}}}
  #  x = t(r %*% eta)
  return(x)
}
