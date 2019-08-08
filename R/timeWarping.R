################################################################################
#
#    Time warping functions
#
################################################################################
#' Srivastava timewarping of one session or two sessions simultanously, both
#' must be measured on the same grid
#'
#' @param DATA1 list containing the Euler angles 3 x N matrix for each trial
#' @param DATA2 list containing the Euler angles 3 x N matrix for each trial
#' @param times vector containing the measurement times
#' @return list/vector with elements the time vector to which the data
#' corresponds.
#' @export
#'
time.warpingSriva <- function( DATA1, DATA2 = NULL, times ){
  if(is.null(DATA2)){
    Yangles  <-  do.call(cbind, lapply( DATA1, function(list) list[2,] ))
    M1 <- M2 <- length(Yangles)
  }else{
    Yangles1  <-  lapply( DATA1, function(list) list[2,] )
    M1 <- length(Yangles1)
    Yangles2  <-  lapply( DATA2, function(list) list[2,] )
    M2 <- M1 + length(Yangles2)
    Yangles   <- cbind( do.call(cbind, Yangles1), do.call(cbind, Yangles2))
  }

  warp     <- time_warping( Yangles*Radian, times, showplot=FALSE )

  if(is.null(DATA2)){
    value <- warp$gam
  }else{
    value <- list( times1 = warp$gam[,1:M1], times2 = warp$gam[,(M1+1):M2] )
  }
  value
}

#' Computes the time warping of curves in SO(3) w.r.t. the total variation loss
#'
#' @param f list geodesicInterpolation object of a curve or if N=NULL a 4xN
#'          matrix containing the data evaluated on a uniform grid of [0,1]
#' @param g list geodesicInterpolation object of a curve.
#' @param N amuount of measurement points of the coarser grid
#' @param factorN2M amount of points between two measurements of the corser grid
#' @return list with elements
#'  \itemize{
#'   \item opt.times Numeric vector containing the optimal grid
#'   \item gamma.index Numeric vector containing the index of the grid points
#'   \item opt.value Numeric optimal value of the minimized functional
#' }
#' @export
#'
timeWarpingSO3 <- function(f, g, N=NULL, factorN2M=2, beta = 2){
  if(is.null(N)){
    N <- length(f$times)
    # get grid size
    M <- N + factorN2M * ( N - 1 )
    # get grid
    fGrid <- seq( 0,1, length.out = N )
    gGrid <- seq( 0,1, length.out = M )
    # get interpolated Data on grids
    Intpol.f <- eval.geodInterpolation( f, fGrid, out = "Quat" )
    Intpol.g <- eval.geodInterpolation( g, gGrid, out = "Quat" )
  }else{
    # get grid size
    M <- N + factorN2M * ( N - 1 )
    # get grid
    fGrid <- seq( 0,1, length.out = N )
    gGrid <- seq( 0,1, length.out = M )

    # get interpolated Data on grids
    Intpol.f <- eval.geodInterpolation( f, fGrid, out = "Quat" )
    Intpol.g <- eval.geodInterpolation( g, gGrid, out = "Quat" )
  }

  # Initialize matrix V.
  V           <- matrix( Inf, N, M )
  V[ 1, ]     <- rep( 0, M )
  V           <- optimV( Intpol.f, Intpol.g, V, gGrid, b = beta, R = fGrid[2] )
  wMin        <- V$wMin
  V           <- V$V

  # calculate gamma
  gamma.index       <- rep( 0, N );
  gamma.index[1]    <- 1
  gamma.index[N]    <- M
  gamma.index[N-1]  <- wMin[N,M] + 1

  for(n in (N - 2):2){
    gamma.index[ n ] <- wMin[ n+1, gamma.index[n+1] ] + 1
  }

  opt.times <- gGrid[ gamma.index ]

  list(
    opt.times   = opt.times,
    gamma.index = gamma.index,
    opt.value   = V[N, which.min(V[N,])]
  )
}

#' Computes the Karcher mean of a sample of curves with respect to the total
#' variation distance of curves in SO(3)
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @param init list geodesicInterpolation object of a curve.
#' @param N amount of measurement points of the coarser grid.
#' @param factorN2M amount of points between two measurements of the corser grid.
#' @param max.iter Numeric amount of iterations.
#' @param err Numeric precision of the iteration.
#' @param show.plot Boolean show plots in each iteration step.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed Karcher mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
KarcherMean <- function( A, init, N = 100, factorN2M = 2, beta = 2, max.iter=5, err = 1e-8, show.plot = FALSE){
  mu.curve <- init
  muGrid <- seq(0,1, length.out=N)

  if(show.plot==TRUE){
    plot(NULL, xlim=c(0,1), ylim=c(-30,70), main = "Raw Data")
    for(i in 1:length(A)){
      matlines( muGrid, t(eval.geodInterpolation(A[[i]], muGrid, out="Euler")) * Radian, col="black" )
    }
    matlines( muGrid, t(eval.geodInterpolation(init, muGrid, out="Euler")) * Radian, col="red", lwd=2 )
  }

  mu    <- eval.geodInterpolation( mu.curve, muGrid, out="Quat")
  gam   <- matrix(NA, length(A), N)
  A.gam <- list()
  count <- 0
  conv  <- 1

  while( count < max.iter & conv == 1 ){
    count <- count + 1
    ## Align each curve to the mean
    for(i in 1:length(A)){
      gam[i,] <- timeWarpingSO3( f = mu.curve, g = A[[i]], N = N, factorN2M = factorN2M, beta = beta)$opt.times
      A.gam[[i]] <- eval.geodInterpolation( A[[i]], gam[i,], out="Quat" )
    }
    ## compute a new mean of the aligned curves
    muNew.curve    <- geodInterpolation( apply(do.call(rbind, A.gam), 2,  function(x)
                                         Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), muGrid)
    muNew          <- eval.geodInterpolation(muNew.curve, muGrid, out="Quat")

    ## break iteration if updated and initial means are close
    if(length.curve(mu, muNew) < err){conv=0}
    print(length.curve(mu, muNew))

    ## plot iteration steps
    if(show.plot==TRUE){
      plot(NULL, xlim=c(0,1), ylim=c(-30,70), main = paste("Aligned Data after",count, "Iterations" ) )
      for(i in 1:length(A)){
        if(i==1){
          matlines(muGrid,t( apply(A.gam[[i]], 2, function(col) Quat2Euler(col)*Radian )  ), col="darkgreen")
        }else{
          matlines(muGrid,t( apply(A.gam[[i]], 2, function(col) Quat2Euler(col)*Radian )  ), col="black")
        }
      }
      matlines(muGrid,t( apply(muNew, 2, function(col) Quat2Euler(col)*Radian )  ), col="red", lwd=2)
      matlines(muGrid,t( apply(mu, 2, function(col) Quat2Euler(col)*Radian )  ), col="blue", lwd=2)
    }
    ## Updated mean curve
    mu <- muNew
    mu.curve <- muNew.curve
  }

  list(mu=mu.curve, gamma = gam)
}

#' Align a set of curves using TV metric
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @param init list geodesicInterpolation object of a curve.
#' @param N amount of measurement points of the coarser grid.
#' @param factorN2M amount of points between two measurements of the corser grid.
#' @param max.iter Numeric amount of iterations.
#' @param err Numeric precision of the iteration.
#' @param show.plot Boolean show plots in each iteration step.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
TimeWarpingTV <- function( A, init, N = 100, factorN2M = 3, beta = 2, max.iter=10,
                           err = 1e-8, show.plot = FALSE){
  ## Initialise variables
  M <- length(A)
  kMean <- NULL
  gamma <- matrix(NA, M, N )
  distMeans <- 1
  times <- seq(0,1, length.out=N)

  ## Compute Karcher mean of the sample
  kMean <- KarcherMean( A, init, N = N, factorN2M = factorN2M, beta = beta,
                        max.iter=max.iter, err = err, show.plot = show.plot )
  print( max(apply(kMean$gam, 1, diff)) )

  gamma <- kMean$gam

#   ## Compute mean warping function
#   kMeanGrid   <- seq(0,1, length.out=N)
#   meanWarp    <- colMeans( kMean$gamma )
#   meanWarpInv <- vapply( kMeanGrid,
#                          function(n) LinInterpolationInverse(n,
#                                                              meanWarp, kMeanGrid),
#                          FUN.VALUE=0 )
#   ## Compute mean curve
#   kMean.geod <- eval.geodInterpolation( kMean$mu, times = meanWarpInv,
#                                         out="Quat")
#   ## Align all samples to centered mean orbit
#   for(i in 1:M){
#     gamma[i,] <- timeWarpingSO3( kMean.geod, A[[i]], N = NULL,
#                                  factorN2M = factorN2M, beta = beta )$opt.times
#   }

  list( mu = kMean$mu, gamma = kMean$gam )
}

#' Linear Interpolation of values of a warping function [0,1] to [0,1]
#'
#' @param gam vector representation of a warping function measured on grid.
#' @param grid vector time points of domain of gam (same length as gam).
#' @param time vector new intermediate time points for evaluating gam.
#' @return vec a vecotr containing the values of gam at the intermediate times
#' @export
#'
LinGam <- function( gam, grid, times){
  dgam <- diff(gam)
  dgrid <- diff(grid)

  vapply(times, function(t)
                  ifelse( max(grid)==t,
                          gam[which(grid==t)],
                          gam[max(which(grid <= t))] + dgam[max(which(grid <= t))] / dgrid[max(which(grid <= t))] * ( t - grid[max(which(grid <= t))] )
                          ),
         FUN.VALUE=1
  )
}