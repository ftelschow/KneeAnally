################################################################################
#
#    Methods for Spatial Alignment of Two Sessions
#
################################################################################
#'  Least Squares esimator to find Rotation R, s.t. A = RB
#'
#' @param A matrix (k x N) containing N vectors in R^k
#' @param B matrix (k x N) containing N vectors in R^k
#' @return list with elements
#'  \itemize{
#'   \item R k x k orthogonal matrix
#'   \item mse mean squared error between A and R%*% B
#'   \item times Vector of measurement times
#' }
#' @export
#'
RotEstim <- function( A, B ){
  k = nrow(A)
  AB <- A%*%t(B)
  SVD <- svd(AB)
  
  if(det(AB) > 0 ){
    S <- S1 <- diag(1,k)
  }else if(abs(det(AB)) < 1e-7){
    if(det(SVD$u %*% t(SVD$v)) > 0){
      S1 <- diag(1,k)
    }else{
      S1 <- diag(c(rep(1,k-1),-1))
    }
  }else{
    S <- S1 <- diag(c(rep(1,k-1),-1))
  }
  # Estimates & Output	
  R <- SVD$u %*% S1 %*% t(SVD$v)

  list( R = R, mse = sum(A^2) + sum(B^2) - 2*sum(SVD$d*diag(S)))
}

#'  Using the length loss to find the elements (P,Q) spatially aligning gamma and eta, s.t. $gamma = P eta Q$
#'
#' @param A matrix (k x N) containing N vectors in R^k
#' @param B matrix (k x N) containing N vectors in R^k
#' @return list with elements
#' @export
#'
RotEstimIntrinsic <- function( gamma, eta ){
  Q = optim( c(0,0,0), hatPillQ, m1=gamma, m2=eta )
  P = optim( c(0,0,0), hatPillP, m1=gamma, m2=eta )
  p = Rot2QuatC(Exp.SO3( P$par )); q = Rot2QuatC(Exp.SO3( Q$par ));
  
  gamma=apply(eta, 2, function(x) QuatMultC(p, QuatMultC(x,q)))
  if(gamma[1,1]<0){
    gamma <- -gamma
  }
  
  list( p = p, q=q, gamma=gamma, pval = P$val, qval=Q$val )
}

#' Spatially align two given sessions A and B of knee data in Euler angles using
#' Srivastavas method for the mean curve
#'
#' @param A list geodesicInterpolation object of a session
#' @param B list geodesicInterpolation object of a session
#' @param times vector common time grid for sessions
#' @return list with elements
#'  \itemize{
#'   \item R k x k orthogonal matrix
#'   \item mse mean squared error between A and R%*% B
#'   \item times Vector of measurement times
#' }
#' @export
#'
alignSessions <- function( A, B, times ){
  # Evaluate mean geodesics
  meanAgeodesic <- eval.geodInterpolation(A$mean, times)*Radian
  meanBgeodesic <- eval.geodInterpolation(B$mean, times)*Radian
  
  warpMean <- time.warpingSriva(DATA1=list(meanAgeodesic), DATA2=list(
                                  meanBgeodesic), times=times)

  meanAwarped <- eval.geodInterpolation(A$mean, times=warpMean$times1, 
                                                out="Quat")
  meanBwarped <- eval.geodInterpolation(B$mean, times=warpMean$times2, 
                                                out="Quat")
  
  R<-RotEstim(meanAwarped, meanBwarped)$R
  
  Aeuler <- lapply(as.list(1:length(A$data)), function(l)
                    eval.geodInterpolation(A$data[[l]], times,
                                                               out="Euler"))
  Beuler <- lapply(as.list(1:length(B$data)), function(l)
    apply( eval.geodInterpolation(B$data[[l]], times,
                           out="Quat"),2, function(x) Quat2Euler(R%*%x )) )
  
  warp <- time.warpingSriva(Aeuler, Beuler, times)

  list(Awarp = warp$times1, Bwarp = warp$times2, Ameanwarp=warpMean$times1,
       Bmeanwarp=warpMean$times2, R=R)
}

#' Computes mean curve of a session and warpings and some descriptive statistics
#'
#' @param A list geodesicInterpolation object of a session
#' @param times vector common time grid for sessions
#' @return list with elements
#'  \itemize{
#'   \item R k x k orthogonal matrix
#'   \item mse mean squared error between A and R%*% B
#'   \item times Vector of measurement times
#' }
#' @export
#'
meanAndWarping <- function( A, times){
  ## Get time warping
  warp <- time.warpingSriva(lapply(A, function(list) 
          eval.geodInterpolation(list, times, out="Euler")), times=times)
  
  ## Get warped data, these are now points measured on the grid "times"
  SessWarped <- lapply(as.list(1:length(A)),
                        function(l) eval.geodInterpolation(A[[l]], warp[,l],
                                                           out="Quat"))
  ## Compute the mean curve  
  mean <- geodInterpolation( apply(do.call(rbind, SessWarped), 2,  function(x)
                    Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), times)
  
  list( data=A, warp=warp, mean=mean, times = times)
}

eval.Align <- function(A,B, alignAB){

  Aeuler <- lapply(as.list(1:length(A$data)), function(l)
                      eval.geodInterpolation(A$data[[l]], alignAB$Awarp[,l],
                                                                  out="Euler"))
  Beuler <- lapply(as.list(1:length(B$data)), function(l)
    apply(eval.geodInterpolation(B$data[[l]], alignAB$Bwarp[,l],
                           out="Quat"), 2, function(x) Quat2Euler(alignAB$R%*%x) ) )
  list(A=Aeuler, B=Beuler)
}





#' Compute aligning rotation
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
alignSessionsSO3 <- function( A, Aalign, B, Balign, warping=TRUE, warpingAB=TRUE, factorN2M = 2, beta = 2, iter = 1, ... ){
  R <- matrix(NA, nrow=4, ncol=4)
  ## Compute mean curves
  mA           <- dim(Aalign$gam)[1]
  mB           <- dim(Balign$gam)[1]
  times       <- Aalign$mu$times
  DATA.quatA  <- DATA.quatB <- list()
  
  for( i  in 1 : mA ){
    if( warping ){
      DATA.quatA[[i]] <- eval.geodInterpolation( A[[i]], times = Aalign$gam[i,], out="Quat")    
    }else{
      DATA.quatA[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat") 
    }
  }
  MeanA <-  geodInterpolation( apply(do.call(rbind, DATA.quatA), 2,  function(x)
    Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), times)
  
  for( i  in 1 : mB ){
    if( warping ){
      DATA.quatB[[i]] <- eval.geodInterpolation( B[[i]], times = Balign$gam[i,], out="Quat")    
    }else{
      DATA.quatB[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat") 
    }
  }
  MeanB <-  geodInterpolation( apply(do.call(rbind, DATA.quatB), 2,  function(x)
    Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), times)
  
  if( !warpingAB ){ iter = 1 }
  
  for( i in 1 : iter ){
    if( i > 1){
      MeanB <- geodInterpolation( apply( R%*%eval.geodInterpolation( MeanB, times = warp$opt.times, out="Quat"),2 , Quat2Euler), times )
    }
    ## time warp mean curves
    if( warpingAB ){
        warp <- timeWarpingSO3( MeanA, MeanB, N=dim(MeanA$data)[2], factorN2M = factorN2M, beta = beta )
        timesB <- warp$opt.times
    }else{
        timesB <- times
    }
    ## align warped mean curves
    MeanA.warped <- eval.geodInterpolation( MeanA, times = times, out="Quat")
    MeanB.warped <- eval.geodInterpolation( MeanB, times = timesB, out="Quat")
    
    R <- RotEstim(MeanA.warped, MeanB.warped)$R
    
#     plot(NULL, xlim=c(0,1), ylim=c(-30,70), ...)
#     matlines( times, t( eval.geodInterpolation( MeanA, times = times, out="Euler") )*Radian, lty=1, col = c("darkblue", "darkred", "darkgreen") )
#       matlines( times, t( eval.geodInterpolation( MeanB, times = warp$opt.times, out="Euler") )*Radian, lty=1, col = "darkgrey" )  
#     matlines( times, t( apply( eval.geodInterpolation( MeanB, times = warp$opt.times, out="Quat"), 2, function(col) Quat2Euler(R%*%col) ) )*Radian, lty=1, col = c("cyan", "orange", "green") )
  }
  list( R = R, warp = timesB )
}


#' Compute aligning rotation
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
alignSessionsSO3intrinsic <- function( A, Aalign, B, Balign, warping=TRUE, warpingAB=TRUE, factorN2M = 2, beta = 2, iter = 1, ... ){
  R <- matrix(NA, nrow=4, ncol=4)
  ## Compute mean curves
  mA           <- dim(Aalign$gam)[1]
  mB           <- dim(Balign$gam)[1]
  times       <- Aalign$mu$times
  DATA.quatA  <- DATA.quatB <- list()
  
  for( i  in 1 : mA ){
    if( warping ){
      DATA.quatA[[i]] <- eval.geodInterpolation( A[[i]], times = Aalign$gam[i,], out="Quat")    
    }else{
      DATA.quatA[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Quat") 
    }
  }
  MeanA <-  geodInterpolation( apply(do.call(rbind, DATA.quatA), 2,  function(x)
    Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), times)
  
  for( i  in 1 : mB ){
    if( warping ){
      DATA.quatB[[i]] <- eval.geodInterpolation( B[[i]], times = Balign$gam[i,], out="Quat")    
    }else{
      DATA.quatB[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Quat") 
    }
  }
  MeanB <-  geodInterpolation( apply(do.call(rbind, DATA.quatB), 2,  function(x)
    Quat2Euler(ProjectiveMean( matrix(x,nrow=4))$mean) ), times)
  
  if( !warpingAB ){ iter = 1 }
  
  for( i in 1 : iter ){
    if( i > 1){
      MeanB$data <- apply( MeanB$data, 2, function(x) Quat2Euler(QuatMultC(hatR$p, QuatMultC(Euler2Quat(x),hatR$q))) )
      MeanB$directions <-  lapply( MeanB$directions, function(l) Vec2screwMC( t(Quat2RotC(hatR$q))%*%ScrewM2VecC(l) ) )
    }
    ## time warp mean curves
    if( warpingAB ){
      warp <- timeWarpingSO3( MeanA, MeanB, N=dim(MeanA$data)[2], factorN2M = factorN2M, beta = beta )
      timesB <- warp$opt.times
    }else{
      timesB <- times
    }
    ## align warped mean curves
    MeanA.warped <- eval.geodInterpolation( MeanA, times = times, out="Quat")
    MeanB.warped <- eval.geodInterpolation( MeanB, times = timesB, out="Quat")
    
    hatR <- RotEstimIntrinsic(MeanA.warped, MeanB.warped)
    
    #     plot(NULL, xlim=c(0,1), ylim=c(-30,70), ...)
    #     matlines( times, t( eval.geodInterpolation( MeanA, times = times, out="Euler") )*Radian, lty=1, col = c("darkblue", "darkred", "darkgreen") )
    #       matlines( times, t( eval.geodInterpolation( MeanB, times = warp$opt.times, out="Euler") )*Radian, lty=1, col = "darkgrey" )  
    #     matlines( times, t( apply( eval.geodInterpolation( MeanB, times = warp$opt.times, out="Quat"), 2, function(col) Quat2Euler(R%*%col) ) )*Radian, lty=1, col = c("cyan", "orange", "green") )
  }
  list( R = R, warp = timesB )
}


#' Compute aligning rotation and temporal registration intrinsically and iterative
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
hatPHIpq <- function( f, g, factorN2M = 4, beta = 0, iter = 10, eps = 1e-6, method = "intrinsic", ... ){
  # Set initial conditions for iteration
  diff = 1; i = 0;
  timesf = f$times;    timesg = g$times;
  hatp   = c(1,0,0,0); hatq = c(1,0,0,0);
  f.val  = eval.geodInterpolation( f, times = timesf, out = "Quat");
  gmod   = g
  
  hatp_old = hatp ; hatq_old = hatq ;
  dist_old = distFTC( f.val, eval.geodInterpolation( g, times = timesf, out = "Quat") ) ;
  
  # plot(NULL, xlim=c(0,1), ylim=c(-90,90))
  # matlines(timesg , t(g$data)*180/pi, lty=1, col=1 )
  # matlines(timesf , t(f$data)*180/pi, lty=1, col=2 )
  
  # Theoretically this loss should decrease in every step. It seems that due to rounding error for small angles this does not hold true numerically. Hence we stop iteratiion also if the loss increases.
  while( i < iter && diff > eps ){
    i = i+1;
    
    if( i>1 ){
      dist_old      <- warp$opt.value ;
      hatp_old      <- hatp ;
      hatq_old      <- hatq ;
      timesg_old    <- timesg ; 
    }

    gmod.val <- eval.geodInterpolation( gmod, times = timesg, out = "Quat")
    
    # Estimate spatial alignment between the curves
    if( method == "intrinsic" ){
      hatR <- RotEstimIntrinsic( f.val, gmod.val )
    }else{
      hatR <- Rot42Quat(RotEstimC(f.val, gmod.val))
    }
    hatp <- QuatMultC( hatR$p, hatp )
    hatq <- QuatMultC( hatR$q, hatq )
    
    # update the rotated version of g
    gmod$data <- apply( g$data, 2, function(x) Quat2Euler(QuatMultC(hatp, QuatMultC(Euler2Quat(x),hatq))) )
    gmod$directions <-  lapply( g$directions, function(l) Vec2screwM( t(Quat2Rot(hatq))%*%screwM2Vec(l) ) )
    
    ## time warp mean curves
    warp <- timeWarpingSO3( f, gmod, N = dim(f$data)[2], factorN2M = factorN2M, beta = beta )
    timesg <- warp$opt.times
    plot( NULL, xlim = c(0,1), ylim = c(-90,90), main = i )
    matlines( timesf, t(eval.geodInterpolation( gmod, times = timesg, out = "Euler"))*180/pi, col="red", lty=1 )
    matlines( timesf , t(f$data)*180/pi, lty=1, col=1 )
    
    ## Check convergence
    diff          <- dist_old - warp$opt.value
    print(warp$opt.value)
    print(diff)
  }
  
  if( diff< 0 ){
    hatp <- hatp_old
    hatq <- hatq_old
    timesg <- timesg_old
    gmod$data <- apply( g$data, 2, function(x) Quat2Euler(QuatMultC(hatp, QuatMultC(Euler2Quat(x),hatq))) )
    gmod$directions <-  lapply( g$directions, function(l) Vec2screwM( t(Quat2Rot(hatq))%*%screwM2Vec(l) ) )
  }
  list( hatp = hatp, hatq = hatq, opt.timesg = timesg, gmod = gmod, optvalue = warp$opt.value )
}