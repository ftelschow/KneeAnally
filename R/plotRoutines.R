################################################################################
#
#     Functions for plotting coordinates of knee data
#
################################################################################
#' Plots the raw data
#'
#' @param A List contains the Data
#' @param V Number volunteer
#' @param S Number session
#' @param PlOT Boolian create a new plot window
#' @param side Number limb 1=left, 2=right
#' @param COL vector contains the colors of the three euler angles
#' @export
#' 
plotRaw <- function(A, V, S, side, PLOT=TRUE, COL=c("darkgrey", "pink","lightblue") ){
  DATA <- A[[V,S]][[side]]
  if(PLOT==TRUE){
        plot(NULL,
             xlim = c(0,1),
             ylim = c(min(DATA[[1]]*Radian)-5, max(DATA[[1]]*Radian)+5),
             main = paste("Volunteer", V, "Session", S, "Side", side),
             xlab = "standardized frames",
             ylab = "Euler angles in degree" 
             )
  }
  for(ang in 1:3){
        lapply(DATA, function(l) lines( (0:(length(l[ang,])-1))/(length(l[ang,])-1) ,l[ang,]*Radian, col=COL[ang]) )
  }
}

#' Plots geodesic interpolated data on a given grid
#'
#' @param A List created by geodInterpolation() 
#' @param V Number volunteer
#' @param S Number session
#' @param side Number limb 1=left, 2=right
#' @param PlOT Boolian create a new plot window
#' @param COL vector contains the colors of the three euler angles
#' @export
#' 
plotGeodInt <- function(A, V, S, side, PLOT=TRUE, COL=c("darkgrey", "pink","lightblue") ){
  DATA <- A[[V,S]][[side]]
  if(PLOT==TRUE){
    plot(NULL,
         xlim = c(0,1),
         ylim = c(min(eval.geodInterpolation(DATA[[1]], DATA[[1]]$times)*Radian)-5,
                  max(eval.geodInterpolation(DATA[[1]], DATA[[1]]$times)*Radian+5)),
         main = paste("Volunteer", V, "Session", S, "Side", side),
         xlab = "standardized frames",
         ylab = "Euler angles in degree" 
    )
  }
  for(ang in 1:3){
    lapply(DATA, function(l) lines( l$times, eval.geodInterpolation(l, l$times)[ang,]*Radian, col=COL[ang]) )
  }
}




#' Plots a Trial and its mean after aligning, you may choose weather warping etc
#' is shown.
#'
#' @param angle Number rotation angle in radian
#' @param coord Char either "x","y" or "z", specifying the coordinate axis of
#' the rotation
#' @return Matrix Rotation matrix around the specidied coordinate axis
#' @export
#' 
plotTrials <- function(A, plot=TRUE, mean=TRUE, warp=TRUE, main=NULL,col=c("darkgrey",
                       "pink","lightblue"), colMean=c("black","darkred",
                       "darkblue"), xlim=c(0,1), ylim=c(-25,70), angles=c(1,2,3)
                       ){
  if(plot==TRUE){
    plot(NULL, xlim=xlim, ylim=ylim, main=main, ylab="Euler angles in degree", xlab="relative time")
  }
  if(warp==TRUE){
    Aeuler <- lapply(as.list(1:length(A$data)), function(l)
                      eval.geodInterpolation(A$data[[l]], A$warp[,l],
                                                            out="Euler"))  
  }else{
    Aeuler <- lapply( A$data, function(list) eval.geodInterpolation( list,
                                                                    A$times  ) )
  }

  Aeuler <- do.call(cbind, lapply( Aeuler, function(list) t(list[angles,]) ))*Radian
  matlines(A$times, Aeuler, col=rep(col, dim(Aeuler)[2]/3), lty=1 )

  if(mean == TRUE){
    meanAeuler <- eval.geodInterpolation( A$mean, A$times  )
    matlines(A$times,
             t(meanAeuler[angles,]*Radian), lty=1, lwd=2, col=colMean)
  }
  
}
#' Plots a Trial and its mean after aligning, you may choose wether warping etc
#' ist shown.
#'
#' @param angle Number rotation angle in radian
#' @param coord Char either "x","y" or "z", specifying the coordinate axis of
#' the rotation
#' @return Matrix Rotation matrix around the specidied coordinate axis
#' @export
#' 
plotTrialsABalign <- function(A,B,AB, plot=TRUE, mean=TRUE, colA=c("darkgrey",
                      "pink","lightblue"), colB=c("brown2", "orange","cyan"),
                      colMeanA="darkred", 
                      colMeanB="black", 
                      xlim=c(0,1), ylim=c(-25,70), angles=c(1,2,3), main=NULL
){
  if(plot==TRUE){
    plot(NULL, xlim=xlim, ylim=ylim, main=main, ylab="Euler angles in degree", xlab="relative time")
  }
  Aeuler <- lapply(as.list(1:length(A$data)), function(l)
                    eval.geodInterpolation(A$data[[l]], AB$Awarp[,l],
                           out="Euler"))
  Beuler <- lapply(as.list(1:length(B$data)), function(l)
                    apply( eval.geodInterpolation(B$data[[l]], AB$Bwarp[,l],
                                  out="Quat"),2, function(x) Quat2Euler(AB$R%*%x
                                                                        )) )
 
  for( ang in angles ){
    AeulerAng <- do.call(rbind, lapply( Aeuler, function(list) t(list[ang,]) ))*Radian
    matlines(A$times, t(AeulerAng), col=colA[ang], lty=1 )
    BeulerAng <- do.call(rbind, lapply( Beuler, function(list) t(list[ang,]) ))*Radian
    matlines(B$times, t(BeulerAng), col=colB[ang], lty=1 )
    
    if(mean == TRUE){
      meanAeuler <- eval.geodInterpolation( A$mean, AB$Ameanwarp  )
      lines(A$times,meanAeuler[ang,]*Radian, lty=1, lwd=2, col=colMeanA)
      
      meanBeuler <- apply( eval.geodInterpolation( B$mean, AB$Bmeanwarp, out="Quat"),2,
                           function(x) Quat2Euler(AB$R%*%x) )
      lines(B$times,meanBeuler[ang,]*Radian, lty=1, lwd=2, col=colMeanB)
    }
  }
}

#' Plots a Trial and its mean after aligning, you may choose wether warping etc
#' ist shown.
#'
#' @param angle Number rotation angle in radian
#' @param coord Char either "x","y" or "z", specifying the coordinate axis of
#' the rotation
#' @return Matrix Rotation matrix around the specidied coordinate axis
#' @export
#' 
plotAB <- function( A, B, Aalign, Balign, V, S, SIDE, angles = c(1,2,3),
                   alignAB = FALSE,
                   warping = FALSE,
                   warpingAB = FALSE,
                   confbands = FALSE,
                   alpha = 0.05,
#                    COL1 = c("darkgrey", "pink","lightblue"), 
#                    COL2 = c("grey20", "purple","cyan"),
                   COL1 = rep("lightblue", 3), 
                   COL2 = rep("pink", 3),
                    COLcb1 = rep("blue", 3), 
                    COLcb2 = rep("red", 3),
                   xlim=c(0,1), ylim=NULL,
                   MAIN= paste( "Side", Limb[1], "Speed Walk" )
){
  A <- A[[ V[1], S[1] ]][[ SIDE[1] ]]$data
  B <- B[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
  Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
  Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
  ALIGN <- NULL
  
  if(warpingAB || alignAB){
    ALIGN <- alignSessionsSO3( A, Aalign, B, Balign, warping=warping, warpingAB = warpingAB , factorN2M = 3, beta = 5, iter=1  )
  }

  N1 <- length( A )
  N2 <- length( B )
  Angles <- c( "x", "y", "z" )
  t  <- seq( 0, 1, length.out=100 )
  dt <- seq( 0, 1, length.out=100 )[2]
  
  DATA1 <- DATA2 <- list()
  
  for( i in 1: N1){
    if( warping ){
      times <- Aalign$gamma[i,]
      times <- sort(c(times, diff( times ) / 2 + times[1:(length(times) -1)] ) )
      DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Euler")
    }else{
      times <- seq(0,1,length.out=199)
      DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Euler")      
    }
  }
  for(i in 1:N2){
    if( warping ){
      times <- Balign$gamma[i,]
      times <- sort(c(times, diff( times ) / 2 + times[1:(length(times) -1)] ) )
      DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Euler")
    }else{
      times <- seq(0,1,length.out=199)
      DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Euler")      
    }
  }
  
  if( alignAB ){
    for( i  in 1 : N2 ){
      if( warpingAB ){
          times <- ALIGN$warp
          times <- sort(c(times, diff( times ) / 2 + times[1:(length(times) -1)] ) )
          if( warping ){
            times <- LinGam( gam = Balign$gamma[i,], grid = t, times = times )
          }else{
            times <- LinGam( gam = t, grid = t, times = times )
          }
          DATA2[[i]] <- apply( ALIGN$R %*% eval.geodInterpolation( B[[i]], times = times, out="Quat"), 2, Quat2Euler )
        }else{
          if( warping ){
            times <- Balign$gamma[i,]
            times <- sort(c(times, diff( times ) / 2 + times[1:(length(times) -1)] ) )
          }else{
            times <- seq(0,1,length.out=199)
          }
          DATA2[[i]] <- apply( ALIGN$R %*% eval.geodInterpolation( B[[i]], times = times, out="Quat"), 2, Quat2Euler )
        }
    }
  }else{
    if( warpingAB ){
      for(i in 1:N2){
        if( warping ){
          times <- ALIGN$warp
          times <- sort(c(times, diff( times ) / 2 + times[1:(length(times) -1)] ) )
          if( warping ){
            times <- LinGam( gam = Balign$gamma[i,], grid = t, times = times )
          }else{
            times <- LinGam( gam = t, grid = t, times = times )
          }
          DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Euler")
        }
      }
    }
    
  }
  ## Compute confidence bands of the data
  warpingABvec <- NULL
  ALIGN2 <- ALIGN
  if( warpingAB ){
    warpingABvec <- ALIGN$warp
    if( !alignAB ){
      ALIGN2 <- NULL
    }
  }
  if( confbands ){
  confA <- confBands( A=A, Aalign = Aalign, alpha = alpha, align = NULL,
                      warping = warping, warpingABvec = NULL,
                      show.plot = FALSE )
  confB <- confBands( A=B, Aalign = Balign, alpha = alpha, align = ALIGN2,
                       warping = warping, warpingABvec = warpingABvec,
                       show.plot = FALSE )
  }
  
  ## Plot data and bands
  for( ang in angles){
    if( is.null(ylim) ){
      y.max <- max( c(DATA1[[1]][ang,], DATA2[[1]][ang,]) ) * Radian
      y.min <- min( c(DATA1[[1]][ang,], DATA2[[1]][ang,]) ) * Radian
    }
    if(ang!=2){
      y.max <- y.max + 4
    }
    plot(NULL,
         xlim = c( 0, 1 ),
         ylim = c( y.min - 7.5, y.max + 7.5 ),
         main = MAIN,
         xlab = "standardized frames",
         ylab = paste( Angles[ang],"angle in degree" ) 
    )
    lapply(DATA1, function(l) lines( seq(0,1,dt/2), l[ang,] * Radian, col=COL1[ang]) )
    lapply(DATA2, function(l) lines(  seq(0,1,dt/2),l[ang,] * Radian, col=COL2[ang]) )
    if( confbands ){
    lines( seq(0,1,length.out=100), confA$up[ang,], col=COLcb1, lwd=2 )
    lines( seq(0,1,length.out=100), confA$lo[ang,], col=COLcb1, lwd=2 )
    lines( seq(0,1,length.out=100), confB$up[ang,], col=COLcb2, lwd=2 )
    lines( seq(0,1,length.out=100), confB$lo[ang,], col=COLcb2, lwd=2 )
    }
    legend( "top", paste("Vol", c(V[1],V[2]), "Ses", c(S[1],S[2]) ), col=c( COL1[ang],COL2[ang] ), lty=rep(1,2) )
  }
  if( any(confB$lo > confA$up) || any(confB$up < confA$lo) ){
    print(V[1])
  }
}

plotMeans <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  align     = FALSE,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  show.plot = FALSE,
  show.plot2 = FALSE,
  Snames    = c("A", "B", "C", "D", "E", "F"),
  xlim      = NULL,
  legend.show    = TRUE,
  ylim      = c(-30,70),
  colA      = c("darksalmon", "darksalmon", "darksalmon"),
  colB      = c( "cadetblue2", "cadetblue2", "cadetblue2"),
  colMeanA  = c("red"),
  colMeanB  = c("blue"), ...
){
  if(is.null(Snames)){
    Snames <- c("A", "B", "C", "D", "E", "F")
  }
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
  
  ## Show plot of the test
  par(mfrow = c(3,1),
      oma = c(4.5,0.2,0,0) + 0.2,
      mar = c(0.5,4,1,1) + 0.2)
  for(ang in 1:3){
    if(ang==1){
      plot(NULL,
           xlim = c(0,100),
           ylim = c(min(apply(mean1, 2, Quat2Euler)[ang,]*Radian)-8, max(apply(mean1, 2, Quat2Euler)[ang,]*Radian)+8),
           main = NULL,
           ylab = ylabVec[ang],
           cex      = 1.5,
           xaxt='n', xlab="",
           cex.axis = 1.5,
           cex.lab  = 1.5
      )
    }
    if(ang==2){
      plot(NULL,
           xlim = c(0,100),
           ylim = c(min(apply(mean1, 2, Quat2Euler)[ang,]*Radian)-8, max(apply(mean1, 2, Quat2Euler)[ang,]*Radian)+8),
           main = NULL,
           ylab = ylabVec[ang],
           xaxt='n', xlab="",
           cex      = 1.5,
           cex.axis =1.5,
           cex.lab  = 1.5
      )
    }
    if(ang==3){
      plot(NULL,
           xlim = c(0,100),
           ylim = c(min(apply(mean1, 2, Quat2Euler)[ang,]*Radian)-8, max(apply(mean1, 2, Quat2Euler)[ang,]*Radian)+8),
           main = NULL,
           ylab = ylabVec[ang],
           xlab = "",
           cex      = 1.5,
           cex.axis = 1.5,
           cex.lab  = 1.5
      )
    }  
    
    if(ang==1){
      legend( "top", c("Session C","Session D"), col=c("darkblue","darkred"), lty=rep(1,1), cex=1.4 )
    }
    lines( seq(0,100,length.out=length(mean1[ang,])) ,apply(mean1, 2, Quat2Euler)[ang,]*Radian, col="darkblue")
    lines( seq(0,100,length.out=length(mean1[ang,])) ,apply(mean2, 2, Quat2Euler)[ang,]*Radian, col="darkred")
    
  }
  title(xlab = "percentage of gait cycle in [%]",
        # ylab = "Euler angles",
        outer = TRUE, line = 3,
        cex.lab      = 2)

}