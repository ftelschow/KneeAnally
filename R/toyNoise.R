#' Create toy noise as mixture of sine and cosine.
#' 
#' @export
toyNoise1 <- function( N, timeGrid=seq(0,1,length.out=100), 
                       sigma = rep(1,length(timeGrid)) ){
  t( (sin(pi/2 * timeGrid) %*% t(rnorm(N,0,1)) + 
     cos(pi/2 * timeGrid) %*% t(rnorm(N,0,1))) * sigma)
}

#' Create toy noise from mixture of Gaussians.
#' 
#' @export
toyNoise2 <- function(N, timeGrid=seq(0,1,length.out=100), 
                      sigma = rep(1,length(timeGrid)), nAnchorPoints = 10){
  
  anchorPoints <- seq(from=0, to=1, length.out=nAnchorPoints)
  f <- sapply(anchorPoints, 
              function(pt) dnorm(timeGrid, mean=pt, 
                                 sd=diff(range(timeGrid))/nAnchorPoints))
  fSqSum <- apply(f^2, 1, sum)
  fNorm <- f / sqrt(fSqSum)
  t((fNorm %*% matrix(rnorm(nAnchorPoints * N), nAnchorPoints, N)) * sigma)
}

#' Create toy noise as OU process.
#' 
#' @export
toyNoise3 <- function( N, timeGrid=seq(0,1,length.out=100),
                       sigma = rep(1,length(timeGrid)), alphaOU=5,
                       sigmaOU=sqrt(10), gamma=NULL ){
  ## Initialize matrix containing as rows a path of the process
  Y  = matrix(NA, length(timeGrid), N)
  
  ## Mean Curve vector
  if(is.null(gamma)){
    gamma <- rep(0,length(timeGrid))
  }
  
  ## Get starting values
  Y[1,] <- rnorm(N, mean = 0, sd = sigmaOU / sqrt( 2*alphaOU )) + gamma[1]
  
  ## Compute the Ornstein-Uhlenbeck trajectory
  for (n in 2:length(timeGrid) ){
    dt <- timeGrid[n]-timeGrid[n-1]
    Y[n,] <- vapply( Y[n-1,], function(y) rnorm(1, mean = (y - gamma[n-1])
                              * exp(-alphaOU * dt) + gamma[n], sd = sigmaOU
                              * sqrt((1 - exp(-2 * alphaOU * dt)) / 2 / alphaOU) )
                     , FUN.VALUE=1)
  }
  t(Y * sigma)
}

#' Create toy noise as non Gaussian noise.
#' 
#' @export
toyNoiseNonGauss <- function( N, timeGrid = seq(0,1,length.out=100), sigma = rep(1,length(timeGrid)) ){
  rcoeff <- cbind( rchisq(N, 1, ncp = 0), rexp(N) )
  
  t( vapply( 1:N, function(l) sqrt(2)/6 * (rcoeff[l,1]-1) * sin( pi*timeGrid ) + 2/3 * ( rcoeff[l,2] - 1 ) * ( timeGrid - 0.5 ), FUN.VALUE=rep(0,length(timeGrid)) ) * sigma )
}