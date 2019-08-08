################################################################################
#
#     Different Residuals of SO(3) valued data
#
################################################################################
#' Residuals of the exponential model: y(t) = g(t) Exp(A(t)) using an estimated
#' g(t).
#'
#' @param DATA List Containing 4 x N matrices of the quaternion representations
#'                  of the curves in the rotation group of each trial
#' @param mean Matrix 4 x N matrix of the quaternion representations
#'                  of an estimated mean curve
#' @return residuum Matrix 3 x N matrix containing the residuals 
#' @export
#' 
Exp.residuals <- function(DATA, mean){
  # Amount of trials in DATA
  M <- length(DATA)
  # Amount of sample times
  N <- dim(mean)[2]
  # Initialize residuum
  residuum <- list()
  
  # Loop over trials
  for(m in 1:M){
    val <- vapply( 1:N, function(n) QuatMultC(QuatInvC(mean[,n]), DATA[[m]][,n] ), FUN.VALUE=rep(0,4)  )
    residuum[[m]] <- apply(val, 2, function(x)  screwM2Vec(LogSO3C(Quat2Rot(x)))      )
  }
  # Return residuum
  residuum
}

#' Compares visually the Residuals of two sessions
#'
#' @param DATA List Containing 4 x N matrices of the quaternion representations
#'                  of the curves in the rotation group of each trial
#' @param mean Matrix 4 x N matrix of the quaternion representations
#'                  of an estimated mean curve
#' @return residuum Matrix 3 x N matrix containing the residuals 
#' @export
#' 
CompResPlot <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  show.plot = FALSE,
  xlim      = c(0,1),
  ylim      = NULL, ...
){
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  N1    <- length( A )
  N2    <- length( B )
  maxN  <- max( N1, N2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= N1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = Aalign$gamma[i,], out="Quat")
      }
      if( i <= N2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = Balign$gamma[i,], out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= N1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= N2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  
  ## Get aligning isometry if alignAB == TRUE
  if( alignAB == TRUE ){
    R <- RotEstim( mean1, mean2 )$R
  }else{
    R <- diag(1,4)
  }
  
  ## Get time warping between sessions if warpingAB == TRUE
  if( warpingAB == TRUE ){
    mean1.geod <- geodInterpolation( apply(mean1, 2, Quat2Euler), times  )
    mean2.geod <- geodInterpolation( R %*% mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat") 
  }else{
    mean2 <- R %*% mean2
  }
  
  ## Get the new time warped and spatially aligned data2
  for( i in 1:N2 ){
    if(warpingAB == TRUE){
      t <- timeAB$opt.times  
    }else{
      t <- times
    }
    if( warping == TRUE ){
      t <- LinGam( gam = Balign$gamma[i,], grid = times, times = t )
    }
    DATA2[[i]] <- R %*% eval.geodInterpolation( B[[i]], times = t, out="Quat")
  }
  
  res1 <- Exp.residuals(DATA1, mean1)
  res2 <- Exp.residuals(DATA2, mean2)
  
  if(is.null(ylim)){
    y <- max(abs(cbind(do.call(cbind,res1), do.call(cbind,res2)))) * Radian * 1.15
    ylim <- c(-y,y)
  }
  
  par(mfrow=c(2,1))
  plot(NULL, xlim=xlim, ylim=ylim, ...)
  lapply(res1, function(list) matlines(times, t(list)*Radian, lty=1))
  plot(NULL, xlim=c(0,1), ylim=ylim, ...)
  lapply(res2, function(list) matlines(times, t(list)*Radian, lty=1))
}