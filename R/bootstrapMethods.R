#' Bootstrap Naive approach (Degras 2009,  Simultaneous confidence bands for
#' nonparametric regression with functional data)
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param M Integer (amount of bootstrap replicates)
#' @param alpha Number (quantile of the bootstrap)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the variance will also be bootstrapped.
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap
#'   \item c.alpha alpha-quantile value of z
#' }
#' @export
bootstrapNaive <- function(A, M = 5000, alpha=0.05, sigma=NULL){
  
  A <- t(A)
  distVec <- rep(NA,M)
  n <- dim(A)[2]
  
  # Note that COLUMNS are realizations of rmultinom!
  counter <- rmultinom(M, size=n, prob=rep(1/n, n))
  
  # Compute means, second moments and standard deviations 
  # of bootstrap realizations.
  # Again COLUMNS correspond to realizations.
  bootMeans <- A %*% counter / n
  if(is.null(sigma)){
    bootSecMoments <- A^2 %*% counter / n
    # We put an abs here to make sure that no NaNs are produced due to machine
    # precision error.
    sigma <- sqrt((n / (n-1)) * abs(bootSecMoments - bootMeans^2)) 
  }
  
  meanA <- rowMeans(A)
  
  distVec <- sqrt(n) * apply( abs((bootMeans - meanA) / sigma), 2, max)
  
  list(z=distVec, c.alpha=quantile(distVec,1-alpha, type=8))
}


#' Multiplier Bootstrap approach 
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param M Integer (amount of bootstrap replicates)
#' @param alpha Number (quantile of the bootstrap)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the empirical variance will be used.
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap
#'   \item c.alpha alpha-quantile value of z
#' }
#' @export
bootstrapMultiplier <- function(A, M = 5000, alpha=0.05, sigma=NULL){
  
  A <- t(A)
  n <- dim(A)[2]
  meanA <- rowMeans(A)
  
  if(is.null(sigma)){
    sigma <- sqrt( rowMeans((A-meanA)^2 ) * n / (n-1))
  } 
  
  normA <- (A - meanA) / sigma  
  
  multiplier <- matrix(rnorm(n * M), n, M)
  distVec <- apply( abs(normA %*% multiplier), 2 ,max) / sqrt(n)
  
  list(z=distVec, c.alpha=quantile(distVec,1-alpha, type=8))
}


#' Bootstrap Multiplier approach bootstrapping the variance.
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param M Integer (amount of bootstrap replicates)
#' @param alpha Number (quantile of the bootstrap)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the variance will also be bootstrapped.
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap
#'   \item c.alpha alpha-quantile value of z
#' }
#' @export
bootstrapMultiplierVar <- function(A, M = 5000, alpha=0.05, sigma=NULL){
  
  A <- t(A)  
  n <- dim(A)[2]
  
  meanA <- rowMeans(A)
  normA <- A - meanA
  
  multiplier <- matrix( rnorm(n*M), n, M)
  bootMeans <- normA %*% multiplier / n
  if(is.null(sigma)){
    boot2ndMoments <- (normA^2) %*% (multiplier^2) / n
    sigma <- sqrt((boot2ndMoments - bootMeans^2) / (n-1) * n)
  }
  
  
  distVec <- sqrt(n) * apply( abs(bootMeans / sigma), 2, max)
  
  list(z=distVec, c.alpha=quantile(distVec,1-alpha, type=8))
}

#' Bootstrap Multiplier2 approach bootstrapping the variance.
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param M Integer (amount of bootstrap replicates)
#' @param alpha Number (quantile of the bootstrap)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the variance will also be bootstrapped.
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap
#'   \item c.alpha alpha-quantile value of z
#' }
#' @export
bootstrapMultiplierVar2 <- function(A, M = 5000, alpha=0.05, sigma=NULL){
  
  A <- t(A)  
  n <- dim(A)[2]
  
  meanA <- rowMeans(A)
  normA <- sqrt(n/(n-1))*(A - meanA)
  
  multiplier <- matrix( rnorm(n*M), n, M)
  bootMeans <- normA %*% multiplier / n
  if(is.null(sigma)){
    boot2ndMoments <- (normA^2) %*% (multiplier^2) / n
    sigma <- sqrt((boot2ndMoments - bootMeans^2) / (n-1) * n)
  }
  
  
  distVec <- sqrt(n) * apply( abs(bootMeans / sigma), 2, max)
  
  list(z=distVec, c.alpha=quantile(distVec,1-alpha, type=8))
}

#' Approximate the tail probabilities using the Gaussian kinematic formula.
#' It approximates the tail probability of the absolute value by twice the 
#' usual tail probability of the process.
#' Lipschitz killing curvature is estimated by simple integral approximation.
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param alpha Number (quantile of the bootstrap)
#' @param timeGrid The locations at which the processes are observed. The 
#'                 default is an equally spaced grid on the unit interval with
#'                 number of points equal to the number of columns of A. (Is 
#'                 used if timeGrid=NULL)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the empirical variance will be used.
#' @return The approximated 1-alpha quantile.
#' @export
gaussianKinematicFormula <- function(A, alpha=0.05, 
                                     timeGrid=NULL, 
                                     sigma=NULL){
  
  if(is.null(timeGrid)){
    timeGrid=seq(0,1,length.out=ncol(A))
  }
  
  A <- t(A)
  
  if(is.null(sigma)){
    sigma <- matrixStats::rowSds(A)
  }
  # Pointwise standard deviation of the derivative.
  dSigma <- matrixStats::rowSds(diff(A / sigma) / diff(timeGrid))
  # Integrate dSigma over the domain.
  L1 <- sum(dSigma * diff(timeGrid))
  
  tailProb <- function(u){
    2*(L1 * exp(-u^2 / 2) / (2 * pi) + (1 - pnorm(u)))  - alpha
  }
  
  list(c.alpha = uniroot(tailProb, interval=c(0,5))$root, L1=L1)
}

#' Approximate the tail probabilities using the Gaussian kinematic formula for 
#' t-fields.
#' It approximates the tail probability of the absolute value by twice the 
#' usual tail probability of the process.
#' Lipschitz killing curvature is estimated by hermte functional.
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param alpha Number (quantile of the bootstrap)
#' @param timeGrid The locations at which the processes are observed. The 
#'                 default is an equally spaced grid on the unit interval with
#'                 number of points equal to the number of columns of A. (Is 
#'                 used if timeGrid=NULL)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the empirical variance will be used.
#' @return The approximated 1-alpha quantile.
#' @export
gaussianKinematicFormulaThermite <- function(A, alpha=0.05, 
                                     timeGrid=NULL, 
                                     sigma=NULL, M=1){
  
  if(is.null(timeGrid)){
    timeGrid=seq(0,1,length.out=ncol(A))
  }
  
  # Degrees of freedom.
  nu <- nrow(A) - 1
  
  A <- t(A)
  
  if(is.null(sigma)){
    sigma <- matrixStats::rowSds(A)
  }
  # Pointwise standard deviation of the derivative.
  dSigma <- matrixStats::rowSds(diff(A / sigma) / diff(timeGrid))
  # Integrate dSigma over the domain.
  L1 <- sum(dSigma * diff(timeGrid))
  
  tailProb <- function(u){
    2*( L1 * (1 + u^2 / nu)^((1 - nu) / 2) / (2 * pi) 
       + (1 - pt(u, df=nu)))  - alpha
  }
  
  list(c.alpha = uniroot(tailProb, interval=c(0,50))$root, L1=L1)
}
#' Approximate the tail probabilities using the Gaussian kinematic formula for 
#' t-fields.
#' It approximates the tail probability of the absolute value by twice the 
#' usual tail probability of the process.
#' Lipschitz killing curvature is estimated by hermte functional.
#' 
#' @param A Matrix (the rows are realizations of the random function)
#' @param alpha Number (quantile of the bootstrap)
#' @param timeGrid The locations at which the processes are observed. The 
#'                 default is an equally spaced grid on the unit interval with
#'                 number of points equal to the number of columns of A. (Is 
#'                 used if timeGrid=NULL)
#' @param sigma The pointwise variance of the process. If var=NULL (which is the
#'            default), the empirical variance will be used.
#' @return The approximated 1-alpha quantile.
#' @export
gaussianKinematicFormulaT <- function(A, alpha=0.05, 
                                      timeGrid=NULL, 
                                      sigma=NULL, M=1,
                                      varMethod = "point"){
  
  if(is.null(timeGrid)){
    timeGrid=seq(0,1,length.out=ncol(A))
  }
  
  # Degrees of freedom.
  nu <- nrow(A) - 1
  
  A <- t(A)
  
  if(is.null(sigma)){
    if( varMethod=="point"){
      sigma <- matrixStats::rowSds(A)
    }else{
      sigma <- as.vector( var.shrink(t(A)) )
    }
  }

  # Pointwise standard deviation of the derivative.
  dSigma <- matrixStats::rowSds(diff(A / sigma) / diff(timeGrid))
  # Integrate dSigma over the domain.
  L1 <- sum(dSigma * diff(timeGrid))
  
  tailProb <- function(u){
    2*( L1 * (1 + u^2 / nu)^((1 - nu) / 2) / (2 * pi) 
        + (1 - pt(u, df=nu)))  - alpha
  }
  
  list(c.alpha = uniroot(tailProb, interval=c(0,50))$root, L1=L1)
}

#' Approximate the tail probabilities using the Gaussian kinematic formula for 
#' Hotelling T2-fields.
#' 
#' @param A List (each list contains the residuals at time t as d x n matrix [d dimension of the vector, n sample amount])
#' @param alpha Number (quantile of the bootstrap)
#' @param timeGrid The locations at which the processes are observed. The default is an equally spaced grid on the unit interval with number of points equal to the number of columns of A. (Is used if timeGrid=NULL)
#' @return The approximated 1-alpha quantile.
#' @export
gaussianKinematicFormulaT2 <- function(A, alpha=0.05, 
                                       timeGrid=NULL,
                                       range=300, nu=NULL){
  
  if( is.null(timeGrid) ){
    timeGrid = seq( 0, 1, length.out = length(A) )
  }

  # Degrees of freedom.  
  if( is.null(nu) ){
    nu <- dim( A[[1]] )[2] - 1    
  }

  
  R <- lapply(A, function(list) t(list) )#- rowMeans(list)))
  
  # Estimate L1
  Q <- lapply( R, function(list) apply(list, 2, function(col) col / sqrt(sum(col^2)) )  )
  
  D <- list()
  for( i in 1:(length(A)-1) ){
    D[[i]] <- Q[[i+1]] - Q[[i]]
  }
  
  meanD <- lapply( D, function(list) sum(sqrt(diag(t(list) %*% list))) / dim(D[[1]])[2] )
  L1 =0
  for(k in 1:length(D)){
    L1 <- L1 + meanD[[k]]
  }
  
  tailProb <- function(u){
    (
      L1 * ( (1 + u / nu)^((1 - nu) / 2) / pi + ( -1 + u * (nu-1) / nu ) * (1 + u / nu)^((1 - nu) / 2) / pi ) 
      +      ( 2*(1 - pt(sqrt(u), df=nu))       + 2* gamma((nu+1)/2) * sqrt(u) * (1 + u / nu)^((1 - nu) / 2) / gamma(nu/2) / sqrt( pi*nu ) )
    )  - alpha
  }
  
  list(c.alpha = uniroot(tailProb, interval=c(0,range))$root, L1=L1)
  
}

