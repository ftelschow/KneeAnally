0###################################################################################################
##
##    Simulations of Exp Model
##
###################################################################################################
#### 
#' Function generating exp model samples
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
GenerateDATA <- function( nSamp,
                          times = seq(0,1,length.out=100),
                          gamma0 = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
                          noise="toyNoise2",
                          sigma1vec=rep(0.005,length(times)), sigma2vec=rep(0.005,length(times)), sigma3vec=rep(0.005,length(times)), corr=FALSE,
                          Q = diag(1,3), P = diag(1,3),
                          show.plot=FALSE
){
  # Get noise function
  noiseFunction <- match.fun(noise)
  
  gamma <- lapply( as.list(1:nSamp), function(l) gamma0  )
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  ## Get random samples of gaussian noise and correlate them
  # Generate Noise
  x <- noiseFunction( nSamp, times, sigma = sigma1vec )
  y <- noiseFunction( nSamp, times, sigma = sigma2vec )
  z <- noiseFunction( nSamp, times, sigma = sigma3vec )
  # Correlation if neccessary
  if( corr == TRUE ){
    y1 <- (y + x) / 2
    z <-  (z + y + x) / sqrt(3)
    y <- y1
    x <- x
    rm(y1)
  }
  # glue them together
  A <- lapply(as.list(1:nSamp), function(n) rbind(x[n,],y[n,],z[n,]) )
  
  # Construct data
  DATA <- lapply(as.list(1:nSamp), function(Trial) apply(t(as.matrix(1:length(times))), 2, function(t) Rot2Quat(P%*%Quat2Rot(gamma[[Trial]][,t])%*%Exp.SO3( Vec2screwM(A[[Trial]][,t]))%*%Q) ))
  
  # plot data in Euler
  # par(mai=rep(0,4))
  if(show.plot){
    plot(NULL, xlim=c(0,1), ylim=c(-90,90))
    lapply(DATA, function(l) matlines(times, t(apply(l,2,Quat2Euler))*Radian) )    
  }
  
  DATA
  
}
#################################################################################
#' Simulate Consistency of Rotation estimate
#'
#' @param M Int number of simulations
#' @param nSamp Int Number of samples from exp model
#' @param noise String Choose between "toyNoise1", "toyNoise2" and "toynoise3"
#' @param sigma1 Numeric standard deviation of the noise process
#' @param sigma2 Numeric standard deviation of the noise process
#' @param sigma3 Numeric standard deviation of the noise process
#' @param sigmaName String choose between "const" and "sin" for the variance function depending on t
#' @param timeGrid Vector Time points at which curves are sampled
#' @param corr Boolean correlate the coordinates of the error process
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#'
#################################################################################
Mean.compSim <- function(
  M               = M,
  nSamp           = nSamp,
  alpha           = 0.05,
  noise           = noise,
  sigma1, sigma2, sigma3,
  sigmaName       = "const",
  timeGrid        = seq(0,1,length.out=100),
  Q               = Euler2Rot( c(12,0,5)/180*pi ),
  P               = Euler2Rot( c(-0.5,13,-9)/180*pi ),
  corr            = FALSE
){
  # counter of covering rate
  rate <- 0
  nT <- length(timeGrid)
  
  # Get noise function
  noiseFunction <- match.fun(noise)
  
  # Different models of the error
  if(sigmaName == "const"){
    sigma1vec <- sigma1*rep(1,100)
    sigma2vec <- sigma2*rep(1,100)
    sigma3vec <- sigma3*rep(1,100)
  }else{
    sigma1vec <- sigma1*(sin(4*pi*timeGrid) + 1.5) / 2
    sigma2vec <- sigma2*(sin(4*pi*timeGrid) + 1.5) / 2
    sigma3vec <- sigma3*(sin(4*pi*timeGrid) + 1.5) / 2
  }
  #Generate mean curve
  gamma0 <- vapply( timeGrid, function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t)) / 180 * pi ), FUN.VALUE=rep(0,4))
  
  distVec <- rep(NA,M)
  
  # Simulate trials
  for(m in 1:M){
    ## Get random samples of gaussian noise and correlate them
    # Generate data
    DATA <- GenerateDATA( nSamp, times = timeGrid, noise=noise, gamma=gamma0,
                          sigma1vec=sigma1vec, sigma2vec=sigma2vec, sigma3vec=sigma3vec, corr=corr,
                          Q = Q, P = P,
                          show.plot=FALSE )
    
    mean <- apply(do.call(rbind, DATA), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    
    ## Get aligning isometry if alignAB == TRUE
    R <- RotEstim( gamma0, mean )$R

    ## Correct mean for aligning isometry
    new.mean <- R %*% mean
    
#     plot(NULL, xlim=c(0,1), ylim=c(-70,70))
#     matlines(timeGrid, t(apply(gamma0, 2, Quat2Euler)*180/pi), col="black" )
#     matlines(timeGrid, t(apply(mean, 2, Quat2Euler)*180/pi), col="red" )
#     matlines(timeGrid, t(apply(new.mean, 2, Quat2Euler)*180/pi), col="blue" )
    
    distVec[m] <- sum(new.mean-gamma0)^2
  } # Loop end: Simulate trials
  
  list( 
          distances = distVec,
          data = data.frame( meanDistance = mean(distVec), nSamp = nSamp, noise = noise, f_j = sigmaName,
          sigma     = paste(sigma1,sigma2,sigma3), Correlation=corr
          )
  )
}

#################################################################################
#' Simulate Covering Rate of confidence bands under Exp Model
#'
#' @param M Int number of simulations
#' @param nSamp Int Number of samples from exp model
#' @param noise String Choose between "toyNoise1", "toyNoise2" and "toynoise3"
#' @param sigma1 Numeric standard deviation of the noise process
#' @param sigma2 Numeric standard deviation of the noise process
#' @param sigma3 Numeric standard deviation of the noise process
#' @param sigmaName String choose between "const" and "sin" for the variance function depending on t
#' @param timeGrid Vector Time points at which curves are sampled
#' @param corr Boolean correlate the coordinates of the error process
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#'
#################################################################################
CovRate.ExpModelSim <- function(
  M               = M,
  nSamp           = nSamp,
  alpha           = 0.05,
  noise           = noise,
  sigma1, sigma2, sigma3,
  sigmaName       = "const",
  timeGrid        = seq(0,1,length.out=100),
  corr            = FALSE
){
  # counter of covering rate
  rate <- 0
  nT <- length(timeGrid)
  
  # Get noise function
  noiseFunction <- match.fun(noise)
  
  # Different models of the error
  if(sigmaName == "const"){
    sigma1vec <- sigma1*rep(1,100)
    sigma2vec <- sigma2*rep(1,100)
    sigma3vec <- sigma3*rep(1,100)
  }else if(sigmaName=="atan"){
    sigma1vec <- sigma1*(2 + 2*atan(5*(timeGrid-0.5) )/pi)
    sigma2vec <- sigma2*(2 + 2*atan(5*(timeGrid-0.5) )/pi)
    sigma3vec <- sigma3*(2 + 2*atan(5*(timeGrid-0.5) )/pi)
  }else{
    sigma1vec <- sigma1*(sin(4*pi*timeGrid) + 1.5) / 2
    sigma2vec <- sigma2*(sin(4*pi*timeGrid) + 1.5) / 2
    sigma3vec <- sigma3*(sin(4*pi*timeGrid) + 1.5) / 2
  }
  #Generate mean curve
  gamma0 <- vapply( timeGrid, function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t)) / 180 * pi ), FUN.VALUE=rep(0,4))
  
  # Simulate trials
  for(k in 1:M){
    ## Get random samples of gaussian noise and correlate them
    # Generate data
    DATA <- GenerateDATA( nSamp, times = timeGrid, noise=noise, gamma=gamma0,
                          sigma1vec=sigma1vec, sigma2vec=sigma2vec, sigma3vec=sigma3vec, corr=corr,
                          Q = diag(1,3), P = diag(1,3),
                          show.plot=FALSE )
    
    mean <- apply(do.call(rbind, DATA), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    Y <- Exp.residuals(DATA, mean)
    
    res <- list()
    for( l in 1:length(timeGrid)){
      res[[l]] <- do.call( cbind, lapply(Y, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y)
    ## Sample Covariance matrices
    T1 <- lapply(res, function(list) nSamp * ginv(var(t(list))) )
    
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res, alpha=alpha, timeGrid=timeGrid, range=400)$c.alpha
    
    
    mu.dist <- vapply( 1:nT, function(j) ScrewM2VecC(LogSO3C( t(Quat2RotC(mean[,j]))%*%Quat2Rot(gamma0[,j]) )), FUN.VALUE = rep(0,3))
    
    dist <- vapply(1:nT, function(u) t(mu.dist[,u])%*%T1[[u]]%*%mu.dist[,u], FUN.VALUE=1)
    
    if( all(dist < t1.alpha) ){
      rate <- rate + 1 
    }
  } # Loop end: Simulate trials
  
  data.frame( covRate = rate/M, nSamp = nSamp, noise = noise, f_j = sigmaName,
              sigma = paste(sigma1,sigma2,sigma3), Correlation=corr
  )
}


####################################################################################
#' Simulate Ttest for two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
TtestSim.ExpModel <- function(
  M               = 1000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise1", "toyNoise1"),
  sigma1          = c(0.01, 0.01),
  sigma2          = c(0.01, 0.01),
  sigma3          = c(0.01, 0.01),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec1 <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec1 <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec1 <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec2 <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec2 <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec2 <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations trials
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp = nSamp1, times = times1, gamma0 = gamma0, noise = noise1, sigma1vec = sigma1vec1, sigma2vec = sigma2vec1, sigma3vec = sigma3vec1, corr = corr1, Q = diag(1,3), P = diag(1,3), show.plot = show.plot )
    
    DATA2 <- GenerateDATA( nSamp = nSamp2, times = times2, gamma0 = eta0, noise = noise2, sigma1vec = sigma1vec2, sigma2vec = sigma2vec2, sigma3vec = sigma3vec2, corr = corr2, Q = Q, P = P, show.plot = show.plot )
    
    #### Estimate the means of the sessions
    ## Compute Ziezold mean on the Sphere
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    
    # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
    # lapply(DATA1, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # lapply(DATA2, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }else{
      warpAB = FALSE
    }
    
    ## Estimate time warping
    if( warpAB == TRUE ){
      if( alignAB == TRUE ){
        R1    <- RotEstimC( mean1, mean2 )
        mean2 <- R1 %*% mean2
      }
      
      mean1.geod <- geodInterpolation( mean1, times1  )
      mean2.geod <- geodInterpolation( mean2, times1  )
      timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times1), factorN2M = factorN2M )
      mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
      if( alignAB == TRUE ){
        R2         <- RotEstimC( mean1, mean2 )
        mean2      <- R2 %*% mean2
      }
    }
    # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
    
    ## Get the new time warped and spatially aligned data2
    if( warpAB == TRUE ){
      DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times1) )
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
    
    # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
    # lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
    # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    # lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
    # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    h <- hotTest(DATA1=DATA1, DATA2=DATA2, times=times1, alpha=alpha, show.plot=show.plot)
    
    if( all(h$stat < h$thresh) == TRUE ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
    p.values[m] <- h$p.value
  } # Loop end: Simulate trials
  
  
  list( p.values=p.values, data = data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB ) )
  
}




#################################################################################### 
#' Simulate Permutation test for two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
PermutationTestSim.ExpModel  <- function(
  Mperm           = 5000,
  M               = 2000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise1", "toyNoise1"),
  sigma1          = c(0.01, 0.01),
  sigma2          = c(0.01, 0.01),
  sigma3          = c(0.01, 0.01),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  PQ              = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec1 <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec1 <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec1 <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec2 <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec2 <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec2 <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(PQ==TRUE){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations trials
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp = nSamp1, times = times1, gamma0 = gamma0, noise = noise1, sigma1vec = sigma1vec1, sigma2vec = sigma2vec1, sigma3vec = sigma3vec1, corr = corr1, Q = diag(1,3), P = diag(1,3), show.plot = show.plot )
    
    DATA2 <- GenerateDATA( nSamp = nSamp2, times = times2, gamma0 = eta0, noise = noise2, sigma1vec = sigma1vec2, sigma2vec = sigma2vec2, sigma3vec = sigma3vec2, corr = corr2, Q = Q, P = P, show.plot = show.plot )
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }else{
      warpAB = FALSE
    }

    #### Compute p-value
    p.values[m] <-  PermTest(DATA1, DATA2, times1, Mperm, alignAB, warpAB=warpAB, factorN2M=factorN2M)

    ### Accept/reject?
    if( p.values[m] >= alpha ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
} 

#################################################################################### 
#' Last Hope Permutation test 2: Simulate two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
PermutationTest2Sim.ExpModel  <- function(
  Mperm           = 5000,
  M               = 2000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise1", "toyNoise1"),
  sigma1          = c(0.01, 0.01),
  sigma2          = c(0.01, 0.01),
  sigma3          = c(0.01, 0.01),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec1 <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec1 <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec1 <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec2 <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec2 <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec2 <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations trials
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp = nSamp1, times = times1, gamma0 = gamma0, noise = noise1, sigma1vec = sigma1vec1, sigma2vec = sigma2vec1, sigma3vec = sigma3vec1, corr = corr1, Q = diag(1,3), P = diag(1,3), show.plot = show.plot )
    
    DATA2 <- GenerateDATA( nSamp = nSamp2, times = times2, gamma0 = eta0, noise = noise2, sigma1vec = sigma1vec2, sigma2vec = sigma2vec2, sigma3vec = sigma3vec2, corr = corr2, Q = Q, P = P, show.plot = show.plot )
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }else{
      warpAB = FALSE
    }
    
    #### Compute p-value
    p.values[m] <- Perm2Test(DATA1=DATA1, DATA2=DATA2, times=times1, Mperm=Mperm, warpAB=warpAB, factorN2M=factorN2M)
    
    #### Accept/reject?
    if( p.values[m] >= alpha ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
}  

#################################################################################### 
#' Last Hope Permutation test 3: Simulate two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
PermutationTest3Sim.ExpModel  <- function(
  Mperm           = 1000,
  M               = 2000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise1", "toyNoise1"),
  sigma1          = c(0.01, 0.01),
  sigma2          = c(0.01, 0.01),
  sigma3          = c(0.01, 0.01),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec1 <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec1 <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec1 <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec2 <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec2 <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec2 <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations trials
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp = nSamp1, times = times1, gamma0 = gamma0, noise = noise1, sigma1vec = sigma1vec1, sigma2vec = sigma2vec1, sigma3vec = sigma3vec1, corr = corr1, Q = diag(1,3), P = diag(1,3), show.plot = show.plot )
    
    DATA2 <- GenerateDATA( nSamp = nSamp2, times = times2, gamma0 = eta0, noise = noise2, sigma1vec = sigma1vec2, sigma2vec = sigma2vec2, sigma3vec = sigma3vec2, corr = corr2, Q = Q, P = P, show.plot = show.plot )
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }
    p.values[m] <- Perm3Test( DATA1, DATA2, times1, Mperm, alignAB, warpAB, factorN2M=factorN2M )$pvalue
    
    if( p.values[m] >= alpha ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data = data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
}  


#################################################################################### 
#' Simulate Confidence Band tests for two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#'
####################################################################################
ConfBandsTestSim.ExpModel <- function(
  M               = 1000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise2", "toyNoise2"),
  sigma1          = c(0.05, 0.05),
  sigma2          = c(0.05, 0.05),
  sigma3          = c(0.05, 0.05),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec1 <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec1 <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec1 <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec2 <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec2 <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec2 <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations trials
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp = nSamp1, times = times1, gamma0 = gamma0, noise = noise1, sigma1vec = sigma1vec1, sigma2vec = sigma2vec1, sigma3vec = sigma3vec1, corr = corr1, Q = diag(1,3), P = diag(1,3), show.plot = show.plot )
    
    DATA2 <- GenerateDATA( nSamp = nSamp2, times = times2, gamma0 = eta0, noise = noise2, sigma1vec = sigma1vec2, sigma2vec = sigma2vec2, sigma3vec = sigma3vec2, corr = corr2, Q = Q, P = P, show.plot = show.plot )
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }else{
      warpAB = FALSE
    }
    
    #### Estimate the means of the sessions
    ## Compute Ziezold mean on the Sphere
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    
    # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
    # lapply(DATA1, function(l) matlines( times1, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # lapply(DATA2, function(l) matlines( times1, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    ## Estimate time warping
    if( warpAB == TRUE ){
      if( alignAB == TRUE ){
        R1    <- RotEstimC( mean1, mean2 )
        mean2 <- R1 %*% mean2
      }
      
      mean1.geod <- geodInterpolation( mean1, times1  )
      mean2.geod <- geodInterpolation( mean2, times1  )
      timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times1), factorN2M = factorN2M )
      mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
      if( alignAB == TRUE ){
        R2         <- RotEstimC( mean1, mean2 )
        mean2      <- R2 %*% mean2
      }
    }
    # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
    
    ## Get the new time warped and spatially aligned data2
    if( warpAB == TRUE ){
      DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times1) )
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
    
#     plot(NULL, xlim=c(0,1), ylim=c(-70,70))
#     lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
#     matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
#     lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
#     matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    #### Constructing confidence bands of the means
    ## residuals of the volunteers
    Y1 <- Exp.residuals(DATA1, mean1)
    Y2 <- Exp.residuals(DATA2, mean2)
    
    res1 <- res2 <- list()
    for( l in 1:length(times1)){
      res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
      res2[[l]] <- do.call( cbind, lapply(Y2, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y1)
    rm(Y2)
    ## Sample Covariance matrices
    T1 <- lapply(res1, function(list) nSamp1 * ginv(var(t(list))) )
    T2 <- lapply(res2, function(list) nSamp2 * ginv(var(t(list))) )
    
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times1, range=400)$c.alpha
    t2.alpha <- gaussianKinematicFormulaT2(res2, alpha=alpha, timeGrid=times1, range=400)$c.alpha

    ## Test whether the confidence regions do intersect for all t
    intersec = TRUE
    j        = 0
    # Loop over time points
    while(  j < length(times1) & intersec == TRUE){
      j        <- j+1
      intersec <- FALSE
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      if( t(testVec)%*%T1[[j]]%*%testVec < t1.alpha ){
        intersec <- TRUE
      }
      ## do they intersect on geodesic between the means?
      testVec <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec)) ))
      if(t(testVec)%*%T2[[j]]%*%testVec < t2.alpha){
        intersec <- TRUE
      }
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      if( t(testVec)%*%T2[[j]]%*%testVec < t2.alpha ){
        intersec <- TRUE
      }
      ## do they intersect on geodesic between the means?
      testVec <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec)) ))
      if(t(testVec)%*%T1[[j]]%*%testVec < t1.alpha){
        intersec <- TRUE
      }
      
      ## Compute the eigenvalues at time t of vol1
      SVDu <- svd(T1[[j]])
      ## columns are the eigenvalues of T1[j]
      u <- t(t(SVDu$u) / sqrt(SVDu$d)*sqrt(t1.alpha))
      ## Initialize counter
      e        = 0
      
      while(e < 5e3 & intersec==FALSE){
        ## Sample uniformly from the sphere
        mu <- rnorm(3)
        mu <- mu / sqrt(sum(mu^2))
        ## random sample of ellipsoid
        ue <- mu[1] * u[,1] + mu[2] * u[,2] + mu[3] * u[,3]
        testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(ue)) ))
        
        if(t(testVec)%*%T2[[j]]%*%testVec < t2.alpha){
          intersec <- TRUE
        }
        e <- e + 1
      }
      ## Compute the eigenvalues at time t of vol2
      SVDu <- svd(T2[[j]])
      ## columns are the eigenvalues of T2[j]
      u <- t(t(SVDu$u) / sqrt(SVDu$d)*sqrt(t2.alpha))
      ## Initialize counter
      e        = 0
      while(e < 5e3 & intersec==FALSE){
        ## Sample uniformly from the sphere
        mu <- rnorm(3)
        mu <- mu / sqrt(sum(mu^2))
        ## random sample of ellipsoid
        ue <- mu[1] * u[,1] + mu[2] * u[,2] + mu[3] * u[,3]
        testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(ue)) ))
        
        if(t(testVec)%*%T1[[j]]%*%testVec < t1.alpha){
          intersec <- TRUE
        }
        e <- e + 1
      }
    }# end of loop over times
    
    if( intersec == TRUE ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
}

ConfBandsTest2Sim.ExpModel <- function(
  M               = 1000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise2", "toyNoise2"),
  sigma1          = c(0.05, 0.05),
  sigma2          = c(0.05, 0.05),
  sigma3          = c(0.05, 0.05),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations 
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp=nSamp1, times=times1, gamma0=gamma0, noise = noise1, sigma1vec=sigma1vec1, sigma2vec=sigma2vec1, sigma3vec=sigma3vec1, corr=corr1, Q = diag(1,3), P = diag(1,3), show.plot=show.plot )
    
    DATA2 <- GenerateDATA( nSamp=nSamp2, times=times2, gamma0=eta0, noise = noise2, sigma1vec=sigma1vec2, sigma2vec=sigma2vec2, sigma3vec=sigma3vec2, corr=corr2, Q = Q, P = P, show.plot=show.plot )
    
    #### Apply Hotelling T2 test
    #### Estimate the means of the sessions
    ## Compute Ziezold mean on the Sphere
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    
    
    ## Get aligning isometry if alignAB == TRUE
    if( alignAB == TRUE ){
      R1 <- RotEstim( mean1, mean2 )$R
      ## Correct mean for aligning isometry
      mean2 <- R1 %*% mean2
    }
    
    #     plot(NULL, xlim=c(0,1), ylim=c(-70,70))
    #     matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    #     matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    ## Estimate time warping
    if(all(times1==times2)==FALSE){
      mean1.geod <- geodInterpolation( mean1, times1  )
      mean2.geod <- geodInterpolation( mean2, times2  )
      timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, factorN2M = 5, b=2 )
      mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
    }
    
    ## Get the new time warped and spatially aligned data2
    if(all(times1==times2)==FALSE){
      DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times2  ) )
      t <- timeAB$opt.times
      ## Get the new time warped and spatially aligned data2
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }else{
      DATA2 <- lapply(DATA2, function(l) l)
    }
    
    mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )

    ## Get aligning isometry if alignAB == TRUE
    if( alignAB == TRUE ){
      R1 <- RotEstim( mean12, mean1 )$R
      R2 <- RotEstim( mean12, mean2 )$R
      ## Correct mean for aligning isometry
      mean1 <- R1 %*% mean1
      mean2 <- R2 %*% mean2      
      ## Correct data for aligning isometry
      DATA1 <- lapply(DATA1, function(l) R1%*%l)
      DATA2 <- lapply(DATA2, function(l) R2%*%l)
    }

#         plot(NULL, xlim=c(0,1), ylim=c(-70,70))
#         matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
#         matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
#         matlines(times1, t(apply(mean12, 2, Quat2Euler)*Radian), col="black" )
#         
#             
#         plot(NULL, xlim=c(0,1), ylim=c(-70,70))
#         lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
#         matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
#         lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
#         matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    #### Constructing confidence bands of the means
    ## residuals of the volunteers
    Y1 <- Exp.residuals(DATA1, mean1)
    Y2 <- Exp.residuals(DATA2, mean2)
    
    res1 <- res2 <- list()
    for( l in 1:length(times1)){
      res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
      res2[[l]] <- do.call( cbind, lapply(Y2, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y1)
    rm(Y2)
    ## Sample Covariance matrices
    T1 <- lapply(res1, function(list) nSamp1 * ginv(var(t(list))) )
    T2 <- lapply(res2, function(list) nSamp2 * ginv(var(t(list))) )
    
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times1, range=400)$c.alpha
    t2.alpha <- gaussianKinematicFormulaT2(res2, alpha=alpha, timeGrid=times1, range=400)$c.alpha
    
    ## Test whether the confidence regions do intersect for all t
    intersec = TRUE
    j        = 0
    # Loop over time points
    while(  j < length(times1) & intersec == TRUE){
      j        <- j+1
      intersec <- FALSE
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      if( t(testVec)%*%T1[[j]]%*%testVec < t1.alpha ){
        intersec <- TRUE
      }
      ## do they intersect on geodesic between the means?
      testVec <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec)) ))
      if(t(testVec)%*%T2[[j]]%*%testVec < t2.alpha){
        intersec <- TRUE
      }
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      if( t(testVec)%*%T2[[j]]%*%testVec < t2.alpha ){
        intersec <- TRUE
      }
      ## do they intersect on geodesic between the means?
      testVec <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec)) ))
      if(t(testVec)%*%T1[[j]]%*%testVec < t1.alpha){
        intersec <- TRUE
      }
      
      ## Compute the eigenvalues at time t of vol1
      SVDu <- svd(T1[[j]])
      ## columns are the eigenvalues of T1[j]
      u <- t(t(SVDu$u) / sqrt(SVDu$d)*sqrt(t1.alpha))
      ## Initialize counter
      e        = 0
      
      while(e < 5e3 & intersec==FALSE){
        ## Sample uniformly from the sphere
        mu <- rnorm(3)
        mu <- mu / sqrt(sum(mu^2))
        ## random sample of ellipsoid
        ue <- mu[1] * u[,1] + mu[2] * u[,2] + mu[3] * u[,3]
        testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(ue)) ))
        
        if(t(testVec)%*%T2[[j]]%*%testVec < t2.alpha){
          intersec <- TRUE
        }
        e <- e + 1
      }
      ## Compute the eigenvalues at time t of vol2
      SVDu <- svd(T2[[j]])
      ## columns are the eigenvalues of T2[j]
      u <- t(t(SVDu$u) / sqrt(SVDu$d)*sqrt(t2.alpha))
      ## Initialize counter
      e        = 0
      while(e < 5e3 & intersec==FALSE){
        ## Sample uniformly from the sphere
        mu <- rnorm(3)
        mu <- mu / sqrt(sum(mu^2))
        ## random sample of ellipsoid
        ue <- mu[1] * u[,1] + mu[2] * u[,2] + mu[3] * u[,3]
        testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(ue)) ))
        
        if(t(testVec)%*%T1[[j]]%*%testVec < t1.alpha){
          intersec <- TRUE
        }
        e <- e + 1
      }
    }# end of loop over times
    
    if( intersec == TRUE ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
}

#' Simulate Difference Band tests for two center curves under Exp Model
#'
#' @param nSamp Int Number of samples from exp model
#' @param times Vector Time points at which curves are sampled
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#'
DiffBandsTest2Sim.ExpModel <- function(
  M               = 1000,
  alpha           = 0.05,
  nSamp           = c(10,10),
  gamma0          = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  eta0            = vapply(seq(0,1,length.out=100), function(t) Euler2QuatC(c(80*t^2-80*t+20, t*sin(4*pi*t^(0.7))*70 + 5, 10*cos(13*pi*t))/ 180 * pi ), FUN.VALUE=rep(0,4)),
  noise           = c("toyNoise2", "toyNoise2"),
  sigma1          = c(0.05, 0.05),
  sigma2          = c(0.05, 0.05),
  sigma3          = c(0.05, 0.05),
  sigmaName       = c("const", "const"),
  corr            = c(FALSE, FALSE),
  times           = cbind(seq(0,1,length.out=100),seq(0,1,length.out=100)),
  alignAB         = FALSE,
  factorN2M       = 2,
  show.plot       = FALSE
){
  ## Parameters for generating samples of first distribution
  nSamp1          <- nSamp[ 1]
  noise1          <- noise[ 1]
  times1          <- times[,1]
  nT1             <- length(times1)
  corr1           <- corr[1]
  
  # Different models of the error
  if(sigmaName[1] == "const"){
    sigma1vec1 <- sigma1[1]*rep(1,nT1)
    sigma2vec1 <- sigma2[1]*rep(1,nT1)
    sigma3vec1 <- sigma3[1]*rep(1,nT1)
  }else if(sigmaName[1]=="atan"){
    sigma1vec <- sigma1*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma2vec <- sigma2*(2 + 2*atan(5*(times1-0.5) )/pi)
    sigma3vec <- sigma3*(2 + 2*atan(5*(times1-0.5) )/pi)
  }else{
    sigma1vec1 <- sigma1[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma2vec1 <- sigma2[1]*(sin(4*pi*times1) + 1.5) / 2
    sigma3vec1 <- sigma3[1]*(sin(4*pi*times1) + 1.5) / 2
  }
  
  ## Parameters for generating samples of second distribution
  nSamp2          <- nSamp[ 2]
  noise2          <- noise[ 2]
  times2          <- times[,2]
  nT2             <- length(times2)
  corr2           <- corr[2]
  
  # Different models of the error
  if(sigmaName[2] == "const"){
    sigma1vec2 <- sigma1[2]*rep(1,nT2)
    sigma2vec2 <- sigma2[2]*rep(1,nT2)
    sigma3vec2 <- sigma3[2]*rep(1,nT2)
  }else if(sigmaName[2]=="atan"){
    sigma1vec <- sigma1*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma2vec <- sigma2*(2 + 2*atan(5*(times2-0.5) )/pi)
    sigma3vec <- sigma3*(2 + 2*atan(5*(times2-0.5) )/pi)
  }else{
    sigma1vec2 <- sigma1[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma2vec2 <- sigma2[2]*(sin(4*pi*times2) + 1.5) / 2
    sigma3vec2 <- sigma3[2]*(sin(4*pi*times2) + 1.5) / 2
  }
  
  ## Scaling factor for Euler representation
  Radian        <- 180 / pi
  
  if(alignAB){
    # Rotational misalignment
    Q <- Euler2Rot(c(12,0,5)/Radian)
    P <- Euler2Rot(c(-0.5,13,-9)/Radian)
  }else{
    Q <- P <- diag(1,3)
  }
  ## Acception and p.values
  accept <- rep(NA, M)
  p.values <-  rep(NA, M)
  
  # Loop of simulations 
  for(m in 1:M){
    #### Get two random samples of the exponential models of the two distributions
    DATA1 <- GenerateDATA( nSamp=nSamp1, times=times1, gamma0=gamma0, noise = noise1, sigma1vec=sigma1vec1, sigma2vec=sigma2vec1, sigma3vec=sigma3vec1, corr=corr1, Q = diag(1,3), P = diag(1,3), show.plot=show.plot )
    
    DATA2 <- GenerateDATA( nSamp=nSamp2, times=times2, gamma0=eta0, noise = noise2, sigma1vec=sigma1vec2, sigma2vec=sigma2vec2, sigma3vec=sigma3vec2, corr=corr2, Q = Q, P = P, show.plot=show.plot )
    
    
    #### Get correction for time warping or not
    if( all(times1==times2)==FALSE ){
      warpAB = TRUE
    }else{
      warpAB = FALSE
    }
    
    #### Estimate the means of the sessions
    ## Compute Ziezold mean on the Sphere
    mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
    
    # plot(NULL, xlim=c(0,1), ylim=c(-70,70))
    # lapply(DATA1, function(l) matlines( times1, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # lapply(DATA2, function(l) matlines( times1, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    # matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    ## Estimate time warping
    if( warpAB == TRUE ){
      if( alignAB == TRUE ){
        R1    <- RotEstimC( mean1, mean2 )
        mean2 <- R1 %*% mean2
      }
      
      mean1.geod <- geodInterpolation( mean1, times1  )
      mean2.geod <- geodInterpolation( mean2, times1  )
      timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times1), factorN2M = factorN2M )
      mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")   
      if( alignAB == TRUE ){
        R2         <- RotEstimC( mean1, mean2 )
        mean2      <- R2 %*% mean2
      }
    }
    # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
    
    ## Get the new time warped and spatially aligned data2
    if( warpAB == TRUE ){
      DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times1) )
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
    
    #     plot(NULL, xlim=c(0,1), ylim=c(-70,70))
    #     lapply(DATA1, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="pink" ) )
    #     matlines(times1, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
    #     lapply(DATA2, function(l) matlines(times1, t(apply(l, 2, Quat2Euler)*Radian), col="lightblue" ) )
    #     matlines(times1, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
    
    testN <- min(length(DATA1), length(DATA2))
    #### Constructing confidence bands of the same mean
    diff.DATA <- lapply( as.list(1:testN), function(l) vapply(1:length(times1), function(t) QuatMultC(QuatInvC(DATA1[[l]][,t]), DATA2[[l]][,t]), FUN.VALUE=rep(1,4)) )
    diff.mean <- vapply(1:length(times1), function(p) QuatMultC(QuatInvC(mean1[,p]), mean2[,p] ), FUN.VALUE=rep(0,4))
    ## residuals of the volunteers
    Y1 <- Exp.residuals(diff.DATA, diff.mean)
    
    res1 <- list()
    for( l in 1:length(times1)){
      res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y1)
    ## Sample Covariance matrices
    T1 <- lapply(res1, function(list) nSamp1 * ginv(var(t(list))) )
    
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times1, range=400)$c.alpha
    
    ## Test whether c(1,0,0,0) always included
    intersec = TRUE
    j        = 0
    # Loop over time points
    while(  j < length(times1) & intersec == TRUE){
      j        <- j+1
      intersec = FALSE
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(diff.mean[,j]))%*%diag(3) ))
      if( t(testVec)%*%T1[[j]]%*%testVec < t1.alpha ){
        intersec <- TRUE
      }
    }# end of loop over times
    
    if( intersec == TRUE ){
      accept[m] <- 1
    }else{
      accept[m] <- 0
    }
  } # Loop end: Simulate trials
  
  
  data.frame( accRate = sum(accept) / M, trials = M, alpha=alpha, nSamp1 = nSamp1, nSamp2 = nSamp2, noise1 = noise1, noise2 = noise2, f_j1 = sigmaName[1], f_j2 = sigmaName[2], sigma1 = paste(sigma1[1],sigma2[1],sigma3[1]), sigma2 = paste(sigma1[2],sigma2[2],sigma3[2]), Correlation1=corr1, Correlation2=corr2, alignAB = alignAB )
  
}