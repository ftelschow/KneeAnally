################################################################################
#
#     Confidence bands and test for knee data
#
################################################################################
#' Confidance regions of mean Euler angle curve
#'
#' @param DATA1 list containing the Euler angles 3 x N matrix for each trial
#' @param DATA2 list containing the Euler angles 3 x N matrix for each trial
#' @param times vector containing the measurement times
#' @return list/vector with elements the time vector to which the data
#' corresponds.
#' @export
#'
ConfidanceBands <- function(A, B, AB, alpha=0.05/3,Plot=TRUE,colA=c("darkblue","darkred", "darkgreen"),
                            colB=c("blue","red", "green")){
  results <- rep(0,3)
  nA <- length(A$data)
  nB <- length(B$data)
  
  Aeuler <- lapply(as.list(1:length(A$data)), function(l)
                eval.geodInterpolation(A$data[[l]], AB$Awarp[,l],
                         out="Euler"))
  Beuler <- lapply(as.list(1:length(B$data)), function(l)
    apply(eval.geodInterpolation(B$data[[l]], AB$Bwarp[,l],
                           out="Quat"),2, function(x) Quat2Euler(AB$R%*%x) ))
  
  if(Plot==TRUE){  plot(NULL,xlim=c(0,1),ylim=c(-20,75))
                   legend( "top", c("y-angle Session","x-angle Session","z-angle Session"),
                            col=c("darkblue","green4","pink"),lty=rep(1,3), cex=0.4, 
                           y.intersp=0.2, bty="n" )
  }

  for(a in 1:3){
    angA <- do.call(rbind, lapply( Aeuler, function(list) t(list[a,]) ))*Radian
    meanA  <- colMeans(angA)
    varA   <- sqrt(colMeans(t((meanA-t(angA))^2) * nA / (nA-1)))
    Ac <- gaussianKinematicFormulaT(angA-meanA, alpha=alpha)$c.alpha

    angB <- do.call(rbind, lapply( Beuler, function(list) t(list[a,]) ))*Radian
    meanB  <- colMeans(angB)
    varB   <- sqrt(colMeans(t((meanB-t(angB))^2) * nB / (nB-1)))
    Bc <- gaussianKinematicFormulaT(angB-meanB, alpha=alpha)$c.alpha
    
    upA  <- meanA + Ac*varA /sqrt(nA)
    lowA <- meanA - Ac*varA /sqrt(nA)
    upB  <- meanB + Bc*varB /sqrt(nB)
    lowB <- meanB - Bc*varB /sqrt(nB)
    if(Plot==TRUE){
    lines(A$times, meanA, col=colA[a])
    matlines(A$times, cbind(upA,lowA), col=colA[a],lty=2)
    lines(B$times, meanB, col=colB[a])    
    matlines(B$times, cbind(upB,lowB), col=colB[a],lty=2)
    }
    if(any(upA<lowB) | any(lowA>upB)){results[a] <- which(upA<lowB | lowA>upB)[1]/length(A$times) }
  }
  results
}


#' Compute confidance bands
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
confBands <- function( A, Aalign = NULL, alpha = 0.05, align = NULL, warping = FALSE, warpingABvec = NULL, show.plot = FALSE ){
  nA <- length( A )
  t  <- seq( 0, 1, length.out=100 )
  ##
  if(warping==TRUE){
    time <- Aalign$gamma
  }else{
    time <- matrix(sort(rep(seq(0, 1, length.out=100), nA)), ncol=100)
  }
  ## Get data in euler angles
  DATA.eulerA <- list()
  if( is.null(align) ){
    for( i  in 1 : nA ){
      if( is.null(warpingABvec) ){
        DATA.eulerA[[i]] <- eval.geodInterpolation( A[[i]], times = time[i,], out="Euler")          
      }else{
        time.new <- LinGam( gam = time[i,], grid = t, times = warpingABvec )
        DATA.eulerA[[i]] <- eval.geodInterpolation( A[[i]], times = time.new, out="Euler")
      }
    }
  }else{
    for( i  in 1 : nA ){
      if( is.null(warpingABvec) ){
        DATA.eulerA[[i]] <- apply( align$R %*% eval.geodInterpolation( A[[i]], times = time[i,], out="Quat"), 2, Quat2Euler )
      }else{
        time.new <- LinGam( gam = time[i,], grid = t, time = warpingABvec )
        DATA.eulerA[[i]] <- apply( align$R %*% eval.geodInterpolation( A[[i]], times = time.new, out="Quat"), 2, Quat2Euler )
      }
    }
  }
  
  Angles <- GKF.up <- GKF.lo <- var <- list()
  for( ang in 1:3 ){
    Angles[[ ang ]] <- do.call( cbind, lapply( DATA.eulerA, function( list ) list[ ang, ]*Radian ) )
  }
  
  mean  <- lapply(Angles, rowMeans )
  for( ang in 1 : 3 ){
    var[[ ang ]]   <- sqrt(rowMeans((Angles[[ ang ]] - mean[[ang]])^2) * nA / (nA-1))
  }
  
  for( ang in 1 : 3 ){
    Alpha  <- gaussianKinematicFormulaT( t( (Angles[[ ang ]] - mean[[ ang ]]) / var[[ang]] ), alpha = alpha / 3 )$c.alpha
    
    GKF.up[[ ang ]] <- mean[[ ang ]] + Alpha*var[[ ang ]] / sqrt(nA)
    GKF.lo[[ ang ]] <- mean[[ ang ]] - Alpha*var[[ ang ]] / sqrt(nA)
  }
  GKF.up <- do.call( rbind, GKF.up )
  GKF.lo <- do.call( rbind, GKF.lo )
  mean   <- do.call( rbind, mean )
  
  if( show.plot == TRUE ){
    plot(NULL, xlim=c(0,1), ylim=c(-30,70) )
    matlines( seq(0,1,length.out=100), Angles[[ 1 ]], col="grey", lty = 1 )
    matlines( seq(0,1,length.out=100), Angles[[ 2 ]], col="lightblue", lty = 1 )
    matlines( seq(0,1,length.out=100), Angles[[ 3 ]], col="lightgreen", lty = 1 )
    matlines( seq(0,1,length.out=100), t(mean), col = c("black", "darkblue", "darkgreen"), lty=1 )
    matlines( seq(0,1,length.out=100), t(GKF.up), col = c("black", "darkblue", "darkgreen"), lty=2 )
    matlines( seq(0,1,length.out=100), t(GKF.lo), col = c("black", "darkblue", "darkgreen"), lty=2 )
  }
  list( mean = mean, up.bound = GKF.up, lo.bound = GKF.lo )
}

#' Compute ellipsoidal confidance regions in the Exp Model
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
ConfBandsTest.expModelHot <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  alignShift = FALSE,
  testM     = 1e4,
  factorN2M = 2,
  rejectPlot= TRUE,
  stancePlot= FALSE,
  stance.y = -25,
  show.plot = FALSE,
  show.plot2 = FALSE,
  Snames    = c("A", "B", "C", "D", "E", "F"),
  xlim      = c(0, 100*max(times) ),
  legend.show    = TRUE,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  warpAB <- warpingAB
  
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
  
  # use the mean shift spatial alignment if neccessary
  if( alignShift==TRUE ){
    a = meanAngleAlign(apply( mean1,2, Quat2Euler ), apply( mean2,2, Quat2Euler ))$shift
    
    DATA1 = lapply( DATA1, function(l) apply( l, 2, function(col) Euler2Quat(Quat2Euler(col)-a) ) )
    mean1 = apply( mean1, 2, function(col) Euler2Quat(Quat2Euler(col)-a) )
  }
  
  #### Constructing confidence bands of the same mean
  ## residuals of the volunteers
  Y1 <- Exp.residuals(DATA1, mean1)
  Y2 <- Exp.residuals(DATA2, mean2)
  
  res1 <- res2 <- list()
  for( l in 1:length(times)){
    res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
    res2[[l]] <- do.call( cbind, lapply(Y2, function(list) list[,l]) )
  }
  ## Gebe speicher wieder frei
  rm(Y1)
  rm(Y2)
  ## Sample Covariance matrices
  T1 <- lapply(res1, function(list) nSamp1 * ginv(var(t(list))) )
  T2 <- lapply(res2, function(list) nSamp2 * ginv(var(t(list))) )
  
  ## Sample Covariance matrices
  #   T1 <- lapply(res1, function(list) (N1-1)*var(t(list)) )
  #   T2 <- lapply(res2, function(list) (N2-1)*var(t(list)) )
  #   m1 <- lapply(res1, function(list) rowMeans(list) )
  #   m2 <- lapply(res2, function(list) rowMeans(list) )
  #   H <- vapply(1:length(res1), function(t) (t(m1[[t]] - m2[[t]]))%*%( ginv(T1[[t]]+T2[[t]])* (N1 + N2 - 2) )%*%(m1[[t]] - m2[[t]])*N1*N2/(N1+N2) , FUN.VALUE=1  )  
  
  
  ## p-value estimated using GKF
  t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times, range=400)$c.alpha
  t2.alpha <- gaussianKinematicFormulaT2(res2, alpha=alpha, timeGrid=times, range=400)$c.alpha
  
  ## Test whether the confidence regions do intersect
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- rep(FALSE, length(times))
  # Loop over time points
  for(j in 1:length(times)){
    intersec = FALSE
    ## Test whether the means are contained in the regions
    testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
    if( t(testVec)%*%T1[[j]]%*%testVec < t1.alpha ){
      intersec <- TRUE
    }
    ## do they intersect on geodesic between the means?
    testVec <- testVec / as.vector(sqrt(t(testVec)%*%T1[[j]]%*%testVec))*sqrt(t1.alpha)
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
    testVec <- testVec / as.vector(sqrt(t(testVec)%*%T2[[j]]%*%testVec))*sqrt(t2.alpha)
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
    
    while(e < testM & intersec==FALSE){
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
    ## Compute the eigenvalues at time t of vol1
    SVDu <- svd(T2[[j]])
    ## columns are the eigenvalues of T2[j]
    u <- t(t(SVDu$u) / sqrt(SVDu$d)*sqrt(t2.alpha))
    ## Initialize counter
    e        = 0
    while(e < testM & intersec==FALSE){
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
    # di = svd()$d
    # if(hypothesis2[j]==FALSE){
    #   hypothesis2[j] <- ifelse(di[3]==0, TRUE, FALSE)      
    # }

    hypothesis[j]  <- intersec 
  }# end of loop over times
  ## Show plot of the test
  if(show.plot==TRUE){
    plot(NULL, xlim = xlim, ylim=ylim, ...)
    lapply(DATA1, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colA, lty=1)  )
    lapply(DATA2, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colB, lty=1)  )
    matlines(100*times, t(apply(mean1, 2, Quat2Euler) * Radian), col=colMeanA, lty=1)
    matlines(100*times, t(apply(mean2, 2, Quat2Euler) * Radian), col=colMeanB, lty=1)
    ## show the stance phase by a black bar.
    if(stancePlot==TRUE){
      mean1Euler <- apply(mean1, 2, Quat2Euler)[2,] * Radian
      vel1       <- diff(mean1Euler)
      mean2Euler <- apply(mean2, 2, Quat2Euler)[2,] * Radian
      vel2       <- diff(mean2Euler)
      
      HC1 <- which.min(mean1Euler[1:40]) - 1
      HC2 <- which.min(mean2Euler[1:40]) - 1
      
      HC  <- round( (HC1 + HC2)/2 )
      TO1 <- which.max(vel1[75:100])+74
      TO2 <- which.max(vel2[75:100])+74
      TO  <- round( (TO1+TO2)/2 )
      lines( x = c(HC,TO), y = c(stance.y,stance.y), lwd=2 )
      points(x = c(HC,TO), y = c(stance.y,stance.y), pch="|", cex=2)
      text( (TO-HC)/2 + HC+1 , stance.y -5, labels="Stance Phase", cex=2.0)
    }
    if(rejectPlot == TRUE){
      abline(v=100*times[which(hypothesis==FALSE)])
      if(legend.show){
        legend( "top", inset=.05, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) ), "Reject"), col=c( colMeanA,colMeanB, "black" ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }else{
      if(legend.show){
        legend( x = 27, y= 69, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) )), col=c( colMeanA,colMeanB ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }
    
  }
  mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  Y1 <- Exp.residuals(DATA1, mean12)
  Y2 <- Exp.residuals(DATA2, mean12)
  
  if( show.plot2 ){
    T1bound2 <- T2bound2 <- T1bound1 <- T2bound1 <- matrix(NA, 3, length(times)) 
    for(j in 1:length(times)){
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      T1bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T1bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j])) ))
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      T2bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T2bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j])) ))
      
    }
    Cv1 <- c("darksalmon", "lightgreen", "lightblue")
    Cv2 <- c( "red", "darkgreen", "blue")
    
    for(i in 1:3){
      Ylim <- max( max(abs(Y1[[1]][i,]))+.4*max(abs(Y1[[1]][i,])), max(abs(Y2[[1]][i,]))+.4*max(abs(Y2[[1]][i,])) )
      
      plot(NULL, xlim = c(0, 100*max(times) ), ylim=c(-Ylim,Ylim), ...)
      lapply(Y1, function(list) lines(100*times, list[i,], col= "blue", lty=1)  )
      lapply(Y2, function(list) lines(100*times, list[i,], col= "green", lty=1)  )
      lines(100*times, T2bound1[i,], col= "black",lwd=2, lty=2)
      lines(100*times, T2bound2[i,], col= "black",lwd=2)
      lines(100*times, T1bound1[i,], col= "red",lwd=2, lty=2)
      lines(100*times, T1bound2[i,], col= "red",lwd=2)
      abline(v=100*times[which(hypothesis==FALSE)], h=0)
    }
    
  }
  ## Output the result of the tests    
  hypothesis
}




#' Compute maximum of two sample hottelling T2 statistic for two samples of curves
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
hotTest <- function(DATA1, DATA2, times=seq(0,1,length.out=100), alpha=0.05, show.plot=FALSE, ylim=c(-70,70),...){
  ## Amount of samples in each data set
  N1 <- length(DATA1)
  N2 <- length(DATA2)
  
  ## Pooled sample Mean
  mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  
  
  #### Constructing T Statistic at eacht time point
  ## residuals of the volunteers
  Y1 <- Exp.residuals(DATA1, mean12)
  Y2 <- Exp.residuals(DATA2, mean12)
  
  res1 <- res2 <- list()
  for( l in 1:length(times)){
    res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
    res2[[l]] <- do.call( cbind, lapply(Y2, function(list) list[,l]) )
  }
  ## Gebe speicher wieder frei
  rm(Y1)
  rm(Y2)
  
  ## Sample Covariance matrices
  #   T1 <- lapply(res1, function(list) (N1-1)*var(t(list)) )
  #   T2 <- lapply(res2, function(list) (N2-1)*var(t(list)) )
  #   m1 <- lapply(res1, function(list) rowMeans(list) )
  #   m2 <- lapply(res2, function(list) rowMeans(list) )
  #   H <- vapply(1:length(res1), function(t) (t(m1[[t]] - m2[[t]]))%*%( ginv(T1[[t]]+T2[[t]])* (N1 + N2 - 2) )%*%(m1[[t]] - m2[[t]])*N1*N2/(N1+N2) , FUN.VALUE=1  )
  
  ## Two Sample Hotelling T Statistic of the data at each time point
  H <- vapply(1:length(res1), function(t) hotelling.stat(t(res1[[t]]), t(res2[[t]]))$st , FUN.VALUE=1  )
  
  ## p-value estimated using GKF
  range= 2000
  Res <- lapply(res1, function(list) t(list) )#- rowMeans(list)))
  
  # Estimate L1
  Q <- lapply( Res, function(list) apply(list, 2, function(col) col / sqrt(sum(col^2)) )  )
  
  D <- list()
  for( i in 1:(length(Res)-1) ){
    D[[i]] <- Q[[i+1]] - Q[[i]]
  }
  
  meanD <- lapply( D, function(list) sum(sqrt(diag(t(list) %*% list))) / dim(D[[1]])[2] )
  L1 =0
  for(k in 1:length(D)){
    L1 <- L1 + meanD[[k]]
  }
  
  # t1.L1 <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times, range=range)$L1
  # t2.L1 <- gaussianKinematicFormulaT2(res2, alpha=alpha, timeGrid=times, range=range)$L1
  
  # t3.L1  <- gaussianKinematicFormulaT2(c(res1,res2), alpha=alpha, timeGrid=times, range=range)$L1
  
  nu = N1 + N2 - 2
  
  # L1 = t3.L1
  #   L1 <- t1.L1#max(c(t1.L1,t2.L1))
  
  tailProb <- function(u){
    ( L1 * ( (1 + u / nu)^((1 - nu) / 2) / pi + ( -1 + u * (nu-1) / nu ) * (1 + u / nu)^((1 - nu) / 2) / pi ) 
      +      ( 2*(1 - pt(sqrt(u), df=nu))       + 2* gamma((nu+1)/2) * sqrt(u) * (1 + u / nu)^((1 - nu) / 2) / gamma(nu/2) / sqrt( pi*nu ) )
    )  - alpha
  }
  
  thresh <- uniroot(tailProb, interval=c(0,range))$root
  p.value <- tailProb(max(H)) + alpha
  ## Approximation using expected Euler characteristic may rise over 1
  if(p.value>1){
    p.value <- 1
  }
  
  ## Test whether the confidence regions do intersect
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- vapply(H, function(c) c<thresh, FUN.VALUE = TRUE )
  
  
  ## Show plot of the test
  if( show.plot ){
    plot(NULL, xlim = c(0, 100*max(times) ), ylim=ylim, ...)
    lapply(DATA1, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col= c("darksalmon", "darksalmon", "darksalmon"), lty=1)  )
    lapply(DATA2, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=c( "cadetblue2", "cadetblue2", "cadetblue2"), lty=1)  )
    matlines(100*times, t(apply(mean1, 2, Quat2Euler) * Radian),  col="red", lty=1)
    matlines(100*times, t(apply(mean2, 2, Quat2Euler) * Radian),  col="blue", lty=1)
    matlines(100*times, t(apply(mean12, 2, Quat2Euler) * Radian), col="black", lty=1)
    abline(v=100*times[which(hypothesis==FALSE)])
  }
  list( stat=H, thresh = thresh, p.value = p.value)
}

#' Compute maximum of two sample hottelling T2 statistic for two samples from our data set with aligning etc.
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
hotTestKneeData.expModel <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  rejectPlot= TRUE,
  testM     = 1e4,
  show.plot = FALSE,
  xlim      = NULL,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  Snames <- c("A", "B", "C", "D", "E", "F")
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
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= N2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
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
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times) )
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
      t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = t )
    }
    DATA2[[i]] <- R %*% eval.geodInterpolation( B[[i]], times = t, out="Quat")
  }
  
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
      
      #       R                 <- RotEstim(mean1, data.aligned[,,t])$R
      #       data.aligned[,,t] <- R %*% data.aligned[,,t]
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
  }
  
  H <- hotTest(DATA1=DATA1, DATA2=DATA2, times=times, alpha=alpha, show.plot=TRUE)
  if(show.plot){
    legend( x = 27, y= 69, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) ), "Reject"), col=c( "red","blue", "black" ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=1.6 )
  }
  
  ## Test whether the confidence regions do intersect
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- H$stat >H$thresh
  
  ## Output the result of the tests    
  list(p.value=H$p.value, accept=hypothesis )
}

#' Permutation test on kneeling data
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
PermutationTestFTdistance <- function(dataA = A, dataB = B, Aalign, Balign,
                                      V, S, SIDE,
                                      times     = times,
                                      align     = TRUE,
                                      alignAB   = TRUE,
                                      warping   = FALSE,
                                      warpingAB = TRUE,
                                      factorN2M = 2,
                                      show.plot = FALSE,
                                      Mperm = 2000){
  
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  warpAB <- warpingAB
  
  PermTest( DATA1=DATA1, DATA2=DATA2, times=times, Mperm, alignAB=alignAB, align=align, warpAB=warpingAB, factorN2M=factorN2M )
}

#' Distance of two samples from the data set
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
SessionwiseFTdistance <- function(dataA = A, dataB = B, Aalign, Balign,
                                  V, S, SIDE,
                                  times     = times,
                                  align     = FALSE,
                                  alignAB   = TRUE,
                                  warping   = FALSE,
                                  warpingAB = TRUE,
                                  show.plot = FALSE
){
  
  Snames <- c("A", "B", "C", "D", "E", "F")
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
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= N2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
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
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times) )
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
      t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = t )
    }
    DATA2[[i]] <- R %*% eval.geodInterpolation( B[[i]], times = t, out="Quat")
  }
  
  
  data.aligned <- array(NA, dim=c(4,length(times),N1+N2))
  for(i in 1:N1){
    data.aligned[,,i] <- DATA1[[i]]
  }
  for(i in 1:N2){
    data.aligned[,,N1+i] <- DATA2[[i]]
  }
  ## Spatial alignment for each trial to its mean
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
      
      #       R                 <- RotEstim(mean1, data.aligned[,,t])$R
      #       data.aligned[,,t] <- R %*% data.aligned[,,t]
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
  }
  
  N = N1+N2
  testSamp1 <- data.aligned[,,1:N1]
  testSamp2 <- data.aligned[,,(N1+1):N]
  
  mean1 <- apply( testSamp1, 2,  function(x) ProjectiveMean( x )$mean )
  mean2 <- apply( testSamp2, 2,  function(x) ProjectiveMean( x )$mean )
  
  if(show.plot){
    plot(NULL, main= paste("V1=",V[1]," V2=",V[2], " Session", S, sep = ""), xlim=c(0,1), ylim=c(-40,70) )
    for(i in 1:N){
      if(i <= N1){
        matlines(times, 180 / pi * t(apply(data.aligned[,,i],2, Quat2Euler)), col="pink" )
      }else{
        matlines(times, 180 / pi * t(apply(data.aligned[,,i],2, Quat2Euler)), col="lightblue" )
      }
    }
    matlines(times, 180 / pi * t(apply(mean1,2, Quat2Euler)), col="darkred" )
    matlines(times, 180 / pi * t(apply(mean2,2, Quat2Euler)), col="darkblue" )
  }
  dist.FT(mean1, mean2)
}


#' Permutation test on kneeling data
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
PermutationTest2FTdistance <- function(dataA = A, dataB = B, Aalign, Balign,
                                       V, S, SIDE,
                                       times     = times,
                                       align     = TRUE,
                                       alignAB   = TRUE,
                                       warping   = FALSE,
                                       warpingAB = TRUE,
                                       factorN2M = 2,
                                       show.plot = FALSE,
                                       Mperm = 2000
){
  
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN      <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  Perm2Test(DATA1, DATA2, times, Mperm, warpAB=warpingAB, factorN2M=factorN2M)
}


#' Full Group action Permutation test on kneeling data
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
PermutationTest3FTdistance <- function(dataA = A, dataB = B, Aalign, Balign,
                                       V, S, SIDE,
                                       times     = times,
                                       align     = FALSE,
                                       alignAB   = TRUE,
                                       warping   = FALSE,
                                       warpingAB = TRUE,
                                       show.plot = FALSE,
                                       Mperm = 2000){
  
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN      <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  #### Estimate the means of the sessions
  ## Compute Ziezold mean on the Sphere
  mean1 <- apply(do.call(rbind, DATA1), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  
  ## Get aligning isometry if alignAB == TRUE
  if( alignAB == TRUE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
  }
  mean1.geod <- geodInterpolation( mean1, times  )
  mean2.geod <- geodInterpolation( mean2, times  )
  
  # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
  # matlines(times, t(apply(mean1, 2, Quat2Euler)*Radian), col="red" )
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="blue" )
  
  ## Get time warping between sessions if warpingAB == TRUE
  if( warpingAB == TRUE ){
    mean1.geod <- geodInterpolation( mean1, times  )
    mean2.geod <- geodInterpolation( mean2, times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod, N=length(times), factorN2M = 10 )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat")
    if( alignAB == TRUE ){
      R2    <- RotEstimC( mean1, mean2 )
      mean2 <- R2 %*% mean2
    }
  }
  # matlines(times, t(apply(R2%*%mean2, 2, Quat2Euler)*Radian), col="black" )
  
  ## Get the new time warped and spatially aligned data2
  if( warpingAB == TRUE ){
    t <- timeAB$opt.times
    ## Get the new time warped data2
    DATA22 <- lapply(B, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
  }
  
  mean2 <- apply(do.call(rbind, DATA22), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  ## Get aligning isometry if alignAB == TRUE
  if( alignAB == TRUE ){
    R1    <- RotEstimC( mean1, mean2 )
    mean2 <- R1 %*% mean2
  }
  
  # matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="black" )
  
  
  data.aligned <- array(NA, dim=c(4,length(times),nSamp1+nSamp2))
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
  
  distVec <- PermTest3Loop( Mperm, data.aligned, Perm, PermC, 10, 2 )
  
  # if(show.plot){
  #   plot(NULL, main= paste("V1=",V[1]," V2=",V[2], " Session", S, sep = ""), xlim=c(0,1), ylim=c(-40,70) )
  #   for(i in 1:N){
  #     if(i <= N1){
  #       matlines(times, 180 / pi * t(apply(data.aligned[,,i],2, Quat2Euler)), col="pink" )
  #     }else{
  #       matlines(times, 180 / pi * t(apply(data.aligned[,,i],2, Quat2Euler)), col="lightblue" )
  #     }
  #   }
  #   matlines(times, 180 / pi * t(apply(mean1,2, Quat2Euler)), col="darkred" )
  #   matlines(times, 180 / pi * t(apply(mean2,2, Quat2Euler)), col="darkblue" )
  # }
  list(z = sort(distVec), dist.mean = distFTC(mean1, mean2), pvalue = length(which(sort(distVec)>distFTC(mean1, mean2)))/Mperm )
}

#' Full Group action Permutation test on kneeling data
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
PermutationTest3FTdistanceTry <- function(dataA = A, dataB = B, Aalign, Balign,
                                          V, S, SIDE,
                                          times     = times,
                                          align     = FALSE,
                                          alignAB   = TRUE,
                                          warping   = FALSE,
                                          warpingAB = TRUE,
                                          factorN2M = 2,
                                          show.plot = FALSE,
                                          Mperm = 2000){
  
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN      <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  Perm3Test(DATA1, DATA2, times, Mperm, alignAB, warpAB=warpingAB, factorN2M=factorN2M)
}






#' Permutation test on kneeling data
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
PermutationTest22FTdistance <- function(dataA = A, dataB = B, Aalign, Balign,
                                        V, S, SIDE,
                                        times     = times,
                                        align     = TRUE,
                                        alignAB   = TRUE,
                                        warping   = FALSE,
                                        warpingAB = TRUE,
                                        factorN2M = 2,
                                        show.plot = FALSE,
                                        Mperm = 2000,
                                        stanceCut = FALSE
){
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- max( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  warpAB <- warpingAB
  
  Perm22Test( DATA1=DATA1, DATA2=DATA2, times=times, Mperm, alignAB=alignAB, align=align, warpAB=warpingAB, factorN2M=factorN2M, stanceCut=stanceCut)
}





#' Compute ellipsoidal confidance regions of differences of Exp Models
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
DiffBandsTest.expModelHot <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  factorN2M = 2,
  rejectPlot= TRUE,
  stancePlot= FALSE,
  stance.y = -25,
  show.plot = FALSE,
  show.plot2 = FALSE,
  Snames    = c("A", "B", "C", "D", "E", "F"),
  xlim      = c(0, 100*max(times) ),
  legend.show    = TRUE,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- min( nSamp1, nSamp2 )
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
#  warpAB <- warpingAB
  
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
#  matlines(times, t(apply(mean2, 2, Quat2Euler)*Radian), col="darkgreen" )
  
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
  
  testN <- min(length(DATA1), length(DATA2))
  #### Constructing confidence bands of the same mean
  diff.DATA <- lapply( as.list(1:testN), function(l) vapply(1:length(times), function(t) QuatMultC(QuatInvC(DATA1[[l]][,t]), DATA2[[l]][,t]), FUN.VALUE=rep(1,4)) )
  diff.mean <- vapply(1:length(times), function(p) QuatMultC(QuatInvC(mean1[,p]), mean2[,p] ), FUN.VALUE=rep(0,4))
  ## residuals of the volunteers
  Y1 <- Exp.residuals(diff.DATA, diff.mean)

  res1 <- list()
  for( l in 1:length(times)){
    res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
  }
  ## Gebe speicher wieder frei
  rm(Y1)
  ## Sample Covariance matrices
  T1 <- lapply(res1, function(list) nSamp1 * ginv(var(t(list))) )

  ## Sample Covariance matrices
  #   T1 <- lapply(res1, function(list) (N1-1)*var(t(list)) )
  #   T2 <- lapply(res2, function(list) (N2-1)*var(t(list)) )
  #   m1 <- lapply(res1, function(list) rowMeans(list) )
  #   m2 <- lapply(res2, function(list) rowMeans(list) )
  #   H <- vapply(1:length(res1), function(t) (t(m1[[t]] - m2[[t]]))%*%( ginv(T1[[t]]+T2[[t]])* (N1 + N2 - 2) )%*%(m1[[t]] - m2[[t]])*N1*N2/(N1+N2) , FUN.VALUE=1  )  
  
  
  ## p-value estimated using GKF
  t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times, range=400)$c.alpha
  
  ## Test whether 0 is inside the confidence regions
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- rep(FALSE, length(times))
  # Loop over time points
  for(j in 1:length(times)){
    intersec = FALSE
    ## Test whether the means are contained in the regions
    testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(diff.mean[,j]))%*%diag(3) ))
    if( t(testVec)%*%T1[[j]]%*%testVec < t1.alpha ){
      hypothesis[j] <- TRUE
    }
  }# end of loop over times
  ## Show plot of the test
  if(show.plot==TRUE){
    plot(NULL, xlim = xlim, ylim=ylim)#, ...)
    lapply(DATA1, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colA, lty=1)  )
    lapply(DATA2, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colB, lty=1)  )
    matlines(100*times, t(apply(mean1, 2, Quat2Euler) * Radian), col=colMeanA, lty=1)
    matlines(100*times, t(apply(mean2, 2, Quat2Euler) * Radian), col=colMeanB, lty=1)
    ## show the stance phase by a black bar.
    if(stancePlot==TRUE){
      mean1Euler <- apply(mean1, 2, Quat2Euler)[2,] * Radian
      vel1       <- diff(mean1Euler)
      mean2Euler <- apply(mean2, 2, Quat2Euler)[2,] * Radian
      vel2       <- diff(mean2Euler)
      
      HC1 <- which.min(mean1Euler[1:40]) - 1
      HC2 <- which.min(mean2Euler[1:40]) - 1
      
      HC  <- round( (HC1 + HC2)/2 )
      TO1 <- which.max(vel1[75:100])+74
      TO2 <- which.max(vel2[75:100])+74
      TO  <- round( (TO1+TO2)/2 )
      lines( x = c(HC,TO), y = c(stance.y,stance.y), lwd=2 )
      points(x = c(HC,TO), y = c(stance.y,stance.y), pch="|", cex=2)
      text( (TO-HC)/2 + HC+1 , stance.y -5, labels="Stance Phase", cex=2.0)
    }
    if(rejectPlot == TRUE){
      abline(v=100*times[which(hypothesis==FALSE)])
      if(legend.show){
        legend( "top", inset=.05, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) ), "Reject"), col=c( colMeanA,colMeanB, "black" ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }else{
      if(legend.show){
        legend( x = 27, y= 69, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) )), col=c( colMeanA,colMeanB ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }
    
  }
  mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  Y1 <- Exp.residuals(DATA1, mean12)
  Y2 <- Exp.residuals(DATA2, mean12)
  
  if( show.plot2 ){
    T1bound2 <- T2bound2 <- T1bound1 <- T2bound1 <- matrix(NA, 3, length(times)) 
    for(j in 1:length(times)){
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      T1bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T1bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j])) ))
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      T2bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T2bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j])) ))
      
    }
    Cv1 <- c("darksalmon", "lightgreen", "lightblue")
    Cv2 <- c( "red", "darkgreen", "blue")
    
    for(i in 1:3){
      Ylim <- max( max(abs(Y1[[1]][i,]))+.4*max(abs(Y1[[1]][i,])), max(abs(Y2[[1]][i,]))+.4*max(abs(Y2[[1]][i,])) )
      
      plot(NULL, xlim = c(0, 100*max(times) ), ylim=c(-Ylim,Ylim), ...)
      lapply(Y1, function(list) lines(100*times, list[i,], col= "blue", lty=1)  )
      lapply(Y2, function(list) lines(100*times, list[i,], col= "green", lty=1)  )
      lines(100*times, T2bound1[i,], col= "black",lwd=2, lty=2)
      lines(100*times, T2bound2[i,], col= "black",lwd=2)
      lines(100*times, T1bound1[i,], col= "red",lwd=2, lty=2)
      lines(100*times, T1bound2[i,], col= "red",lwd=2)
      abline(v=100*times[which(hypothesis==FALSE)], h=0)
    }
    
  }
  ## Output the result of the tests    
  hypothesis
}


#' Compute ellipsoidal confidance regions of differences of Exp Models
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
DiffBandsTestBoots.expModelHot <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  Mboots     = 500,
  factorN2M = 2,
  rejectPlot= TRUE,
  stancePlot= FALSE,
  stance.y = -25,
  show.plot = FALSE,
  show.plot2 = FALSE,
  Snames    = c("A", "B", "C", "D", "E", "F"),
  xlim      = c(0, 100*max(times) ),
  legend.show    = TRUE,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- max( nSamp1, nSamp2 )
  Tt <- length(times)
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
#  warpAB <- warpingAB
  
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
  
  ## Get the new time warped data2
  if( warpAB == TRUE ){
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    if( alignAB == TRUE ){
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }else{
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  DATA1array <- array(0, dim=c(4,Tt,nSamp1))
  DATA2array <- array(0, dim=c(4,Tt,nSamp2))
  meanarray  <- array(0, dim=c(4,Tt,Mboots))
  T1         <- array(0,dim=c(3,3,Tt))
  
  for(u in 1:maxN){
    if(u<=nSamp1){
    DATA1array[,,u] <- DATA1[[u]]
    }
    if(u<=nSamp2){
      DATA2array[,,u] <- DATA2[[u]]
    }
  }
  
  testN <- min(nSamp1, nSamp2)  
  for(tau in 1:Mboots){
    sampleVec1 <- sample.int( nSamp1, size = testN, replace = TRUE, prob = NULL )
    sampleVec2 <- sample.int( nSamp2, size = testN, replace = TRUE, prob = NULL )
    boot1 <- DATA1array[,,sampleVec1]
    boot2 <- DATA2array[,,sampleVec2]
    
    bmean1 <- vapply(1:Tt, function(z) ProjectiveMeanC( boot1[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
    bmean2 <-  vapply(1:Tt, function(z) ProjectiveMeanC( boot2[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
    
    # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
    # apply(boot1, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # apply(boot2, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times, t(apply(bmean1, 2, Quat2Euler)*Radian), col="red", lwd=2 )
    # matlines(times, t(apply(bmean2, 2, Quat2Euler)*Radian), col="blue", lwd=2 )
    
    Rboot      <- RotEstimC( bmean1, bmean2 )
    bmean2     <- Rboot %*% mean2
    # matlines(times, t(apply(bmean2, 2, Quat2Euler)*Radian), col="green", lwd=2 )
    
    boot2      <- array(unlist( mapply(function(u) Rboot %*% boot2[,,u], 1:testN, SIMPLIFY = FALSE)), dim = c(4, Tt, testN))
    
    # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
    # apply(boot1, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # apply(boot2, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times, t(apply(bmean1, 2, Quat2Euler)*Radian), col="red", lwd=2 )
    # matlines(times, t(apply(bmean2, 2, Quat2Euler)*Radian), col="blue", lwd=2 )
    
    diff.bDATA <- lapply( as.list(1:testN), function(l) vapply(1:Tt, function(t) QuatMultC(QuatInvC(boot1[,t,l]), boot2[,t,l]), FUN.VALUE=rep(1,4)) )
    diff.bmean <- vapply(1:length(times), function(p) QuatMultC(QuatInvC(bmean1[,p]), bmean2[,p] ), FUN.VALUE=rep(0,4))
    ## residuals of the volunteers
    Y1 <- Exp.residuals(diff.bDATA, diff.bmean)
    res1 <- list()
    for( l in 1:length(times)){
      res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y1)
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times, range=400)$c.alpha
    ## Sample Covariance matrices
    T1 <- T1 + array(unlist(lapply(res1, function(list) nSamp1 * ginv(t1.alpha*var(t(list))) )), dim=c(3,3,Tt))
    meanarray[,,tau] <- diff.bmean
  }
  T1 <- T1/Mboots
  dmean <- vapply(1:Tt, function(z) ProjectiveMeanC( meanarray[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
  
  ## Test whether 0 is inside the confidence regions
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- rep(FALSE, Tt)
  # Loop over time points
  for(j in 1:Tt){
    intersec = FALSE
    ## Test whether the means are contained in the regions
    testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(dmean[,j]))%*%diag(3) ))
    if( t(testVec)%*%T1[,,j]%*%testVec < 1 ){
      hypothesis[j] <- TRUE
    }
  }# end of loop over times
  ## Show plot of the test
  if(show.plot==TRUE){
    plot(NULL, xlim = xlim, ylim=ylim)#, ...)
    lapply(DATA1, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colA, lty=1)  )
    lapply(DATA2, function(list) matlines(100*times, t(apply(R2%*%R1%*%list, 2, Quat2Euler) * Radian), col=colB, lty=1)  )
    matlines(100*times, t(apply(mean1, 2, Quat2Euler) * Radian), col=colMeanA, lty=1)
    matlines(100*times, t(apply(R2%*%R1%*%mean2, 2, Quat2Euler) * Radian), col=colMeanB, lty=1)
    ## show the stance phase by a black bar.
    if(stancePlot==TRUE){
      mean1Euler <- apply(mean1, 2, Quat2Euler)[2,] * Radian
      vel1       <- diff(mean1Euler)
      mean2Euler <- apply(mean2, 2, Quat2Euler)[2,] * Radian
      vel2       <- diff(mean2Euler)
      
      HC1 <- which.min(mean1Euler[1:40]) - 1
      HC2 <- which.min(mean2Euler[1:40]) - 1
      
      HC  <- round( (HC1 + HC2)/2 )
      TO1 <- which.max(vel1[75:100])+74
      TO2 <- which.max(vel2[75:100])+74
      TO  <- round( (TO1+TO2)/2 )
      lines( x = c(HC,TO), y = c(stance.y,stance.y), lwd=2 )
      points(x = c(HC,TO), y = c(stance.y,stance.y), pch="|", cex=2)
      text( (TO-HC)/2 + HC+1 , stance.y -5, labels="Stance Phase", cex=2.0)
    }
    if(rejectPlot == TRUE){
      abline(v=100*times[which(hypothesis==FALSE)])
      if(legend.show){
        legend( "top", inset=.05, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) ), "Reject"), col=c( colMeanA,colMeanB, "black" ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }else{
      if(legend.show){
        legend( x = 27, y= 69, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) )), col=c( colMeanA,colMeanB ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }
    
  }
  mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  Y1 <- Exp.residuals(DATA1, mean12)
  Y2 <- Exp.residuals(DATA2, mean12)
  
  if( show.plot2 ){
    T1bound2 <- T2bound2 <- T1bound1 <- T2bound1 <- matrix(NA, 3, length(times)) 
    for(j in 1:length(times)){
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      T1bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T1bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j])) ))
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      T2bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T2bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j])) ))
      
    }
    Cv1 <- c("darksalmon", "lightgreen", "lightblue")
    Cv2 <- c( "red", "darkgreen", "blue")
    
    for(i in 1:3){
      Ylim <- max( max(abs(Y1[[1]][i,]))+.4*max(abs(Y1[[1]][i,])), max(abs(Y2[[1]][i,]))+.4*max(abs(Y2[[1]][i,])) )
      
      plot(NULL, xlim = c(0, 100*max(times) ), ylim=c(-Ylim,Ylim), ...)
      lapply(Y1, function(list) lines(100*times, list[i,], col= "blue", lty=1)  )
      lapply(Y2, function(list) lines(100*times, list[i,], col= "green", lty=1)  )
      lines(100*times, T2bound1[i,], col= "black",lwd=2, lty=2)
      lines(100*times, T2bound2[i,], col= "black",lwd=2)
      lines(100*times, T1bound1[i,], col= "red",lwd=2, lty=2)
      lines(100*times, T1bound2[i,], col= "red",lwd=2)
      abline(v=100*times[which(hypothesis==FALSE)], h=0)
    }
    
  }
  ## Output the result of the tests    
  hypothesis
}


#' Compute ellipsoidal confidance regions of differences of Exp Models
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
DiffBandsTestPerm.expModelHot <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  alpha     = 0.05,
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  Mboots     = 500,
  factorN2M = 2,
  rejectPlot= TRUE,
  stancePlot= FALSE,
  stance.y = -25,
  show.plot = FALSE,
  show.plot2 = FALSE,
  Snames    = c("A", "B", "C", "D", "E", "F"),
  xlim      = c(0, 100*max(times) ),
  legend.show    = TRUE,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  Snames <- c("A", "B", "C", "D", "E", "F")
  ## Load the correct data of volunteer and session
  A <- dataA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
  B <- dataB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ## For the plots
  Ses <- c("A", "B", "C", "D", "E", "F")
  ## Basic informations of the data    
  nSamp1    <- length( A )
  nSamp2    <- length( B )
  maxN  <- max( nSamp1, nSamp2 )
  Tt <- length(times)
  # arrays for the data
  DATA1 <- DATA2 <- list()
  
  ## Evaluate the geodesic interpolated data either on a constant grid or
  ## on the time registered grid
  if(warping == TRUE){
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        t <- LinGam( gam = Aalign$gamma[i,], grid = seq(0,1,length.out=length(Aalign$gamma[i,])), times = times )
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = t, out="Quat")
      }
      if( i <= nSamp2 ){
        t <- LinGam( gam = Balign$gamma[i,], grid = seq(0,1,length.out=length(Balign$gamma[i,])), times = times )
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = t, out="Quat")
      }
    }
  }else{
    for( i in 1:maxN ){
      if( i <= nSamp1 ){
        DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat")  
      }
      if( i <= nSamp2 ){
        DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat")  
      }
    }
  }
  
  warpAB <- warpingAB
  
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
  
  ## Get the new time warped data2
  if( warpAB == TRUE ){
    DATA2.geod <- lapply( DATA2, function(l) geodInterpolation(l, times) )
    t <- timeAB$opt.times
    ## Get the new time warped data2
    if( alignAB == TRUE ){
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }else{
      DATA2 <- lapply(DATA2.geod, function(l) eval.geodInterpolation( l, times = t, out="Quat") )
    }
    mean2 <- apply(do.call(rbind, DATA2), 2,  function(x) ProjectiveMeanC( matrix(x,nrow=4), MaxIt=100, err=1e-8) )
  }
  
  DATA1array <- array(0, dim=c(4,Tt,nSamp1+nSamp2))
  meanarray  <- array(0, dim=c(4,Tt,Mboots))
  T1         <- array(0, dim=c(3,3,Tt))
  
  for(u in 1:(nSamp1+nSamp2) ){
    if(u<=nSamp1){
      DATA1array[,,u] <- DATA1[[u]]
    }
    if(u>nSamp1){
      DATA1array[,,u] <- DATA2[[u-nSamp1]]
    }
  }
  
  testN <- min(nSamp1, nSamp2)  
  for(tau in 1:Mboots){
    sampleVec1 <- sample.int( nSamp1+nSamp1, size = testN, replace = TRUE, prob = NULL )
    sampleVec2 <- sample.int( nSamp2+nSamp1, size = testN, replace = TRUE, prob = NULL )
    
    sVec11 <- sampleVec1[which( sampleVec1<=nSamp1 )]
    sVec12 <- sampleVec1[which( sampleVec1>nSamp1 )]
    
    if(length(sVec11)>0){
      boot11 <- array( 1, dim=c(4,Tt, length(sVec11)) )
      boot11[,,1:length(sVec11)] <- DATA1array[,,sVec11]
    }else{
      boot11 <- NULL
    }
    
    if(length(sVec12)>0){
      boot12 <- array( NA, dim=c(4,Tt, length(sVec12)) )
      boot12[,,1:length(sVec12)] <- DATA1array[,,sVec12]
    }else{
      boot12 <- NULL
    }
    
    if( is.null(boot11)==TRUE ){
      boot1 <- boot12
    }else if(is.null(boot12)==TRUE){
      boot1 <- boot11
    }else{
      if(dim(boot11)[3]>1){
        bmean11 <- vapply(1:Tt, function(z) ProjectiveMeanC( boot11[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
      }else{
      bmean11 <- boot11[,,1]
    }
      if(dim(boot12)[3]>1){
        bmean12 <-  vapply(1:Tt, function(z) ProjectiveMeanC( boot12[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
      }else{
        bmean12 <- boot12[,,1]
      }
      
      Rboot1      <- RotEstimC( bmean11, bmean12 )
      boot12      <- array(unlist( mapply(function(u) Rboot1 %*% boot12[,,u], 1:dim(boot12)[3], SIMPLIFY = FALSE)), dim = c(4, Tt, dim(boot12)[3]))
      boot1 <- abind(boot11, boot12,rev.along=1)
    }
    
    sVec21 <- sampleVec2[which( sampleVec2<=nSamp1 )]
    sVec22 <- sampleVec2[which( sampleVec2>nSamp1 )]
    
    if( length(sVec21)>0 ){
      boot21 <- array( 1, dim=c(4,Tt, length(sVec21)) )
      boot21[,,1:length(sVec21)] <- DATA1array[,,sVec21]
    }else{
      boot21 <- NULL
    }
    
    if( length(sVec22)>0 ){
      boot22 <- array( 1, dim=c(4,Tt, length(sVec22)) )
      boot22[,,1:length(sVec22)] <- DATA1array[,,sVec22]
    }else{
      boot22 <- NULL
    }
    
    if( is.null(boot21)==TRUE ){
      boot2 <- boot22
    }else if(is.null(boot22)==TRUE){
      boot2 <- boot21
    }else{
      if(dim(boot21)[3]>1){
        bmean21 <- vapply(1:Tt, function(z) ProjectiveMeanC( boot21[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
      }else{
        bmean21 <- boot21[,,1]
      }
      if(dim(boot22)[3]>1){
        bmean22 <-  vapply(1:Tt, function(z) ProjectiveMeanC( boot22[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
      }else{
        bmean22 <- boot22[,,1]
      }
      
      Rboot2      <- RotEstimC( bmean21, bmean22 )
      boot22      <- array(unlist( mapply(function(u) Rboot2 %*% boot22[,,u], 1:dim(boot22)[3], SIMPLIFY = FALSE)), dim = c(4, Tt, dim(boot22)[3]))
      boot2 <- abind(boot21, boot22, rev.along=1)
    }
    
    bmean1 <- vapply(1:Tt, function(z) ProjectiveMeanC( boot1[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
    bmean2 <-  vapply(1:Tt, function(z) ProjectiveMeanC( boot2[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
    
    # plot(NULL, xlim=c(0,1), ylim=c(-20,60))
    # apply(boot1, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="pink" ) )
    # apply(boot2, 3, function(l) matlines( times, t(apply(l,2, Quat2Euler)*180/pi ), col="lightblue" ) )
    # matlines(times, t(apply(bmean1, 2, Quat2Euler)*Radian), col="red", lwd=2 )
    # matlines(times, t(apply(bmean2, 2, Quat2Euler)*Radian), col="blue", lwd=2 )
    
    # Rboot      <- RotEstimC( bmean1, bmean2 )
    # bmean2     <- Rboot %*% mean2
    # # matlines(times, t(apply(bmean2, 2, Quat2Euler)*Radian), col="green", lwd=2 )
    # 
    # boot2      <- array(unlist( mapply(function(u) Rboot %*% boot2[,,u], 1:testN, SIMPLIFY = FALSE)), dim = c(4, Tt, testN))
    
    diff.bDATA <- lapply( as.list(1:testN), function(l) vapply(1:Tt, function(t) QuatMultC(QuatInvC(boot1[,t,l]), boot2[,t,l]), FUN.VALUE=rep(1,4)) )
    diff.bmean <- vapply(1:length(times), function(p) QuatMultC(QuatInvC(bmean1[,p]), bmean2[,p] ), FUN.VALUE=rep(0,4))
    ## residuals of the volunteers
    Y1 <- Exp.residuals(diff.bDATA, diff.bmean)
    res1 <- list()
    for( l in 1:length(times)){
      res1[[l]] <- do.call( cbind, lapply(Y1, function(list) list[,l]) )
    }
    ## Gebe speicher wieder frei
    rm(Y1)
    ## p-value estimated using GKF
    t1.alpha <- gaussianKinematicFormulaT2(res1, alpha=alpha, timeGrid=times, range=400)$c.alpha
    ## Sample Covariance matrices
    T1 <- T1 + array(unlist(lapply(res1, function(list) nSamp1 * ginv(t1.alpha*var(t(list))) )), dim=c(3,3,Tt))
    meanarray[,,tau] <- diff.bmean
  }
  
  T1 <- T1/Mboots
  dmean <- vapply(1:Tt, function(z) ProjectiveMeanC( meanarray[,z,], MaxIt=100, err=1e-8), FUN.VALUE=rep(1,4))
  
  ## Test whether 0 is inside the confidence regions
  # vector collecting where the null hypothesis is TRUE or FALSE
  hypothesis <- rep(FALSE, Tt)
  # Loop over time points
  for(j in 1:Tt){
    intersec = FALSE
    ## Test whether the means are contained in the regions
    testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(dmean[,j]))%*%diag(3) ))
    if( t(testVec)%*%T1[,,j]%*%testVec < 1 ){
      hypothesis[j] <- TRUE
    }
  }# end of loop over times
  ## Show plot of the test
  if(show.plot==TRUE){
    plot(NULL, xlim = xlim, ylim=ylim, ...)
    lapply(DATA1, function(list) matlines(100*times, t(apply(list, 2, Quat2Euler) * Radian), col=colA, lty=1)  )
    lapply(DATA2, function(list) matlines(100*times, t(apply(R2%*%R1%*%list, 2, Quat2Euler) * Radian), col=colB, lty=1)  )
    matlines(100*times, t(apply(mean1, 2, Quat2Euler) * Radian), col=colMeanA, lty=1)
    matlines(100*times, t(apply(R2%*%R1%*%mean2, 2, Quat2Euler) * Radian), col=colMeanB, lty=1)
    ## show the stance phase by a black bar.
    if(stancePlot==TRUE){
      mean1Euler <- apply(mean1, 2, Quat2Euler)[2,] * Radian
      vel1       <- diff(mean1Euler)
      mean2Euler <- apply(mean2, 2, Quat2Euler)[2,] * Radian
      vel2       <- diff(mean2Euler)
      
      HC1 <- which.min(mean1Euler[1:40]) - 1
      HC2 <- which.min(mean2Euler[1:40]) - 1
      
      HC  <- round( (HC1 + HC2)/2 )
      TO1 <- which.max(vel1[75:100])+74
      TO2 <- which.max(vel2[75:100])+74
      TO  <- round( (TO1+TO2)/2 )
      lines( x = c(HC,TO), y = c(stance.y,stance.y), lwd=2 )
      points(x = c(HC,TO), y = c(stance.y,stance.y), pch="|", cex=2)
      text( (TO-HC)/2 + HC+1 , stance.y -5, labels="Stance Phase", cex=2.0)
    }
    if(rejectPlot == TRUE){
      abline(v=100*times[which(hypothesis==FALSE)])
      if(legend.show){
        legend( "top", inset=.05, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) ), "Reject"), col=c( colMeanA,colMeanB, "black" ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }else{
      if(legend.show){
        legend( x = 27, y= 69, c(paste("Vol", c(V[1],V[2]), "Ses", c(Snames[S[1]],Snames[S[2]]) )), col=c( colMeanA,colMeanB ), lty=rep(1,2), bty="o", box.col="white", bg="white", cex=2.2 )}
    }
    
  }
  mean12 <- apply(do.call(rbind, c(DATA1,DATA2) ), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
  Y1 <- Exp.residuals(DATA1, mean12)
  Y2 <- Exp.residuals(DATA2, mean12)
  
  if( show.plot2 ){
    T1bound2 <- T2bound2 <- T1bound1 <- T2bound1 <- matrix(NA, 3, length(times)) 
    for(j in 1:length(times)){
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean1[,j]))%*%Quat2Rot(mean2[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T1[[j]]%*%testVec)*sqrt(t1.alpha)
      T1bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T1bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean1[,j])) ))
      ## Test whether the means are contained in the regions
      testVec <- screwM2Vec(Log.SO3( t(Quat2Rot(mean2[,j]))%*%Quat2Rot(mean1[,j]) ))
      ## do they intersect on geodesic between the means?
      testVec1 <- testVec / sqrt(t(testVec)%*%T2[[j]]%*%testVec)*sqrt(t2.alpha)
      T2bound1[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j]))%*%Exp.SO3(Vec2screwM(testVec1)) ))
      T2bound2[,j] <- screwM2Vec(Log.SO3( t(Quat2Rot(mean12[,j]))%*%(Quat2Rot(mean2[,j])) ))
      
    }
    Cv1 <- c("darksalmon", "lightgreen", "lightblue")
    Cv2 <- c( "red", "darkgreen", "blue")
    
    for(i in 1:3){
      Ylim <- max( max(abs(Y1[[1]][i,]))+.4*max(abs(Y1[[1]][i,])), max(abs(Y2[[1]][i,]))+.4*max(abs(Y2[[1]][i,])) )
      
      plot(NULL, xlim = c(0, 100*max(times) ), ylim=c(-Ylim,Ylim), ...)
      lapply(Y1, function(list) lines(100*times, list[i,], col= "blue", lty=1)  )
      lapply(Y2, function(list) lines(100*times, list[i,], col= "green", lty=1)  )
      lines(100*times, T2bound1[i,], col= "black",lwd=2, lty=2)
      lines(100*times, T2bound2[i,], col= "black",lwd=2)
      lines(100*times, T1bound1[i,], col= "red",lwd=2, lty=2)
      lines(100*times, T1bound2[i,], col= "red",lwd=2)
      abline(v=100*times[which(hypothesis==FALSE)], h=0)
    }
    
  }
  ## Output the result of the tests    
  hypothesis
}