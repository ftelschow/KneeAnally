#' Testing two Data sets using ILL Permutation test
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
PermTest <- function(DATA1, DATA2, times, Mperm, alignAB, align=FALSE, warpAB, factorN2M=2){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # lapply(DATA1, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
  # lapply(DATA2, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  if( warpAB == TRUE ){
    if( alignAB == TRUE ){
        R1    <- RotEstimC( mean1, mean2 )
        mean2 <- R1 %*% mean2
    }
    
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = factorN2M )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
    if( alignAB == TRUE ){
        R2         <- RotEstimC( mean1, mean2 )
        mean2      <- R2 %*% mean2
    }
  }
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
  
  ## Get the new time warped and spatially aligned data2
  if( warpAB == TRUE ){
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    if( alignAB == TRUE ){
      DATA2 <- lapply(DATA2.geod, function(l) R2%*%R1%*% eval.geodInterpolation( l, times = t, out="Quat") )
    }else{
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  if( alignAB == TRUE && warpAB == FALSE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
    DATA2 <- lapply(DATA2, function(l) R1 %*% l )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  if( alignAB ==TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
  }
  
  data.aligned <- array( NA, dim=c(4, length(times), nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }
  
  ### Start ISA
  if( align == TRUE ){
    ## Spatial alignment for each trial to its mean
    data.aligned <- array(NA, dim=c(4,length(times),N1+N2))
    for(i in 1:N1){
      data.aligned[,,i] <- DATA1[[i]]
    }
    for(i in 1:N2){
      data.aligned[,,N1+i] <- DATA2[[i]]
    }
    for(t in 1:N1){
      pqoptim <- optim(rep(0,3), minDistL2, m1=data.aligned[,,t], m2=mean1)
      P <- Exp.SO3( Vec2screwM(pqoptim$par) )
      data.aligned[,,t] <- apply(data.aligned[,,t], 2, function(col) Rot2Quat(P%*%Quat2RotC(col)%*%t(P)) )
      
            # R                 <- RotEstim(mean1, data.aligned[,,t])$R
            # data.aligned[,,t] <- R %*% data.aligned[,,t]
    }
    
    for(t in 1:N2){
      pqoptim <- optim(rep(0,3), minDistL2, m1=data.aligned[,,t+N1], m2=mean2)
      P <- Exp.SO3( Vec2screwM(pqoptim$par) )
      data.aligned[,,t+N1] <- apply(data.aligned[,,t+N1], 2, function(col) Rot2Quat(P%*%Quat2RotC(col)%*%t(P)) )
      
      #       R                     <- RotEstim(mean2, data.aligned[,,t+N1])$R
      #       data.aligned[,,N1+t]  <- R %*% data.aligned[,,t+N1]
    }        
    
    for(i in 1:N1){
      for(i in 1:N1){
        DATA1[[i]] <- data.aligned[,,i]
      }
      for(i in 1:N2){
        DATA2[[i]] <- data.aligned[,,N1+i]
      }    
    }
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    if( alignAB ==TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
    }
  }
  ###### End ISA
  
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # apply(data.aligned[,,1:nSamp1], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  # apply(data.aligned[,,(nSamp1+1):( nSamp1+nSamp2)], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="blue" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="red" )
  
  ## Initialize vector for the distances
  distVec <- rep(NA, Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  Perm  <- vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=FALSE, prob = rep(1/N,N))), 1:nSamp1 )
  PermC <- vapply(1:Mperm, function(t) sort(x[is.na(pmatch(x,Perm[,t]))]), 1:nSamp2 )
  
  p2r <- PermTestLoop(Mperm, data.aligned, Perm, PermC)  
  distVec <- sort(p2r)
  
  length(which(distVec>distFTC(mean1, mean2))) / Mperm
}


#' Testing two Data sets using ILL Permutation test with spatial correction for gait similarities
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
Perm2Test <- function(DATA1, DATA2, times, Mperm, warpAB, factorN2M=2){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # lapply(DATA1, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
  # lapply(DATA2, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  if( warpAB == TRUE ){
    if( alignAB == TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
    }
    
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = factorN2M )

    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )

    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  R1    <- RotEstimC( mean1, mean2 )
  mean2 <- R1 %*% mean2
  
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
  
  data.aligned <- array( NA, dim=c(4, length(times), nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }

  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # apply(data.aligned[,,1:nSamp1], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  # apply(data.aligned[,,(nSamp1+1):( nSamp1+nSamp2)], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="blue" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="red" )
  
  ## Initialize vector for the distances
  distVec <- rep(NA, Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  Perm  <- vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=FALSE, prob = rep(1/N,N))), 1:nSamp1 )
  PermC <- vapply(1:Mperm, function(t) sort(x[is.na(pmatch(x,Perm[,t]))]), 1:nSamp2 )
  
  p2r <- PermTest2Loop(Mperm, data.aligned, Perm, PermC)  
  distVec <- sort(p2r)
  
  length(which(distVec>distFTC(mean1, mean2))) / Mperm
}

#'
Boot2Test <- function(DATA1, DATA2, times, Mperm, warpAB, factorN2M=2){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # lapply(DATA1, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
  # lapply(DATA2, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  if( warpAB == TRUE ){
    if( alignAB == TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
    }
    
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = factorN2M )
    
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  R1    <- RotEstimC( mean1, mean2 )
  mean2 <- R1 %*% mean2
  
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
  
  data.aligned <- array( NA, dim=c(4, length(times), nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # apply(data.aligned[,,1:nSamp1], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  # apply(data.aligned[,,(nSamp1+1):( nSamp1+nSamp2)], 3, function(l) matlines(times, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="blue" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="red" )
  
  ## Initialize vector for the distances
  distVec <- rep(NA, Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  boot1 = vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=TRUE, prob = rep(1/N,N))), 1:nSamp1 )
  boot2 = vapply(1:Mperm, function(t) sort(sample(1:N,nSamp2,replace=TRUE, prob = rep(1/N,N))), 1:nSamp2 )
  
  
  p2r <- PermTest2Loop(Mperm, data.aligned, boot1, boot2)  
  distVec <- sort(p2r)
  
  length(which(distVec>distFTC(mean1, mean2))) / Mperm
}

#' Testing two Data sets using ILL Permutation test with full correction for gait similarities
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
Perm3Test <- function(DATA1, DATA2, times, Mperm, alignAB, warpAB, factorN2M){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  
  ## Get aligning isometry if alignAB == TRUE
  R1    <- RotEstimC( mean1, mean2 )
  mean2 <- R1 %*% mean2
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  mean1.geod <- geodInterpolation( mean1, times  )
  mean2.geod <- geodInterpolation( mean2, times  )
  timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = factorN2M )
  mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
  R2    <- RotEstimC( mean1, mean2 )
  mean2 <- R2 %*% mean2
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="green" )
  
  # ## Get the new time warped and spatially aligned data2
  # DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times  ) )
  # t <- timeAB$opt.times
  # ## Get the new time warped data2
  # DATA22 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
  # mean2 <- apply(do.call(rbind, DATA22), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  #   if( alignAB == TRUE ){
  #     R1    <- RotEstimC( mean1, mean2 )
  #     mean2 <- R1 %*% mean2
  #   }
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  data.aligned <- array( NA, dim=c(4, length(times), nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }
  
  ## Initialize vector for the distances
  distVec <- rep(NA,Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  Perm  <- vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=FALSE, prob = rep(1/N,N))), 1:nSamp1 )
  PermC <- vapply(1:Mperm, function(m) sort(x[is.na(pmatch(x,Perm[,m]))]), 1:nSamp2 )
  
  p2r <- PermTest3Loop(Mperm, data.aligned, Perm, PermC, N2Mfac=factorN2M, b=2)  
  
  list( pvalue = length(which( sort(p2r)>distFTC(mean1, mean2) )) / Mperm, distVec=sort(p2r), distMeans = distFTC(mean1, mean2) )
}



#' Testing two Data sets using ILL Permutation test with full correction for gait similarities
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
Perm3Test.old <- function(DATA1, DATA2, times, Mperm, alignAB, warpAB){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  
  ## Get aligning isometry if alignAB == TRUE
  if( alignAB == TRUE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
  }
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  if( warpAB == TRUE ){
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = 5 )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
    if( alignAB == TRUE ){
      R2    <- RotEstimC( mean1, mean2 )
      mean2 <- R2 %*% mean2
    }
  }
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="green" )
  
  ## Get the new time warped and spatially aligned data2
  if( warpAB == TRUE ){
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times  ) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    DATA22 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    mean2 <- apply(do.call(rbind, DATA22), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    if( alignAB == TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
    }
  }
  
  # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  # lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  data.aligned <- array( NA, dim=c(4, length(times), nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }
  
  ## Initialize vector for the distances
  distVec <- rep(NA,Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  Perm  <- vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=FALSE, prob = rep(1/N,N))), 1:nSamp1 )
  PermC <- vapply(1:Mperm, function(m) sort(x[is.na(pmatch(x,Perm[,m]))]), 1:nSamp2 )
  
  p2r <- PermTest3Loop(Mperm, data.aligned, Perm, PermC, 2, 2)  
  
  list( pvalue=length(which( sort(p2r)>distFTC(mean1, mean2) )) / Mperm, distVec = sort(p2r), distMeans = distFTC(mean1, mean2) )
}






#' Testing two Data sets using ILL Permutation test with spatial correction for gait similarities
#'
#' @param A list containing geodesicInterpolation object of a curve.
#' @return list with elements
#'  \itemize{
#'   \item mu List Interpolation object of the computed center of orbit Karcher
#'         mean.
#'   \item gamma Matrix containing the warping functions to the Karcher mean.
#' }
#' @export
#'
Perm22Test <- function(DATA1, DATA2, times, Mperm, alignAB, align=FALSE, warpAB, factorN2M=2, stanceCut=FALSE){
  nSamp1 <- length(DATA1)
  nSamp2 <- length(DATA2)
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
  # lapply(DATA1, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
  # lapply(DATA2, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Estimate time warping
  if( warpAB == TRUE ){
    if( alignAB == TRUE ){
      R1    <- RotEstimC( mean1, mean2 )
      mean2 <- R1 %*% mean2
    }
    
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = factorN2M )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
    if( alignAB == TRUE ){
      R2         <- RotEstimC( mean1, mean2 )
      mean2      <- R2 %*% mean2
    }
  }
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
  
  ## Get the new time warped and spatially aligned data2
  if( warpAB == TRUE ){
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    if( alignAB == TRUE ){
      DATA2 <- lapply(DATA2.geod, function(l) R2%*%R1%*% eval.geodInterpolation( l, times = t, out="Quat") )
    }else{
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  if( alignAB == TRUE && warpAB == FALSE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
    DATA2 <- lapply(DATA2, function(l) R1 %*% l )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  if( alignAB ==TRUE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
  }
  
  
  if(stanceCut==TRUE){
    mean1Euler <- apply(mean1, 2, Quat2Euler)[2,] * Radian
    vel1       <- diff(mean1Euler)
    mean2Euler <- apply(mean2, 2, Quat2Euler)[2,] * Radian
    vel2       <- diff(mean2Euler)
    
    HC1 <- which.min(mean1Euler[1:40]) - 1
    HC2 <- which.min(mean2Euler[1:40]) - 1
    
    HC <- round( (HC1 + HC2)/2 )
    thresh <- round(HC*0.1)
    HC <- HC - thresh
    TO1 <- which.max(vel1[75:100])+74
    TO2 <- which.max(vel2[75:100])+74
    # TO1 <- which.min(mean1Euler[75:100])+74
    # TO2 <- which.min(mean1Euler[75:100])+74
    TO  <- round( (TO1+TO2)/2 )#+thresh
    
    times1 <- times[1:HC]#, times[TO:length(times)])
    mean1 <- mean1[,1:HC]#,TO:length(times))]

    DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t[1:HC], out="Quat") )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    
    DATA1 <- lapply(DATA1, function(l) l[,1:HC])
    
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
    DATA2 <- lapply( DATA2, function(l) R1%*%l )
  }
  if(stanceCut==TRUE){
    tt <- length(times1)
  }else{
    tt <- length(times)
  }
  
  data.aligned <- array( NA, dim=c(4, tt, nSamp1+nSamp2) )
  for(i in 1:nSamp1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:nSamp2){
    data.aligned[,,nSamp1+i] <- DATA2[[i]]
  }
  
  plot(NULL, xlim=c(0,1), ylim=c(-70,70))
  apply(data.aligned[,,1:nSamp1], 3, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
  apply(data.aligned[,,(nSamp1+1):( nSamp1+nSamp2)], 3, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
  matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="blue" )
  matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="red" )
  
  ## Initialize vector for the distances
  distVec <- rep(NA, Mperm)
  N = nSamp1 + nSamp2
  x <- 1:N
  
  Perm  <- vapply(1:Mperm, function(t) sort(sample(1:N,nSamp1,replace=FALSE, prob = rep(1/N,N))), 1:nSamp1 )
  PermC <- vapply(1:Mperm, function(t) sort(x[is.na(pmatch(x,Perm[,t]))]), 1:nSamp2 )

  distVec <- PermTestLoop(Mperm, data.aligned, Perm, PermC)
  # distVec <- p2r
  # distVec <-  sort(PermTestLoop(Mperm, data.aligned[,(HC+1):dim(mean1)[2],], Perm, PermC))
  length(which(distVec>distFTC(mean1, mean2)))/ Mperm#+
  # length(which(distVec>distFTC( mean1[,(HC+1):dim(mean1)[2]], mean2[,(HC+1):dim(mean2)[2]]))) / Mperm
}
