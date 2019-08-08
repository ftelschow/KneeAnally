################################################################################
##
##    Auswertungsdateien
##
################################################################################

TestAB <- function(
                 A, B, Aalign, Balign,
                 V, S, SIDE,
                 angles    = c(1,2,3),
                 alpha     = 0.05,
                 alignAB   = TRUE,
                 warping   = FALSE,
                 warpingAB = TRUE,
                 confbands = TRUE,
                 MEAN      = "Euler",
                 show.plot = TRUE,
                 COL1      = rep("lightblue", 3), 
                 COL2      = rep("pink", 3),
                 COLcb1    = rep("blue", 3), 
                 COLcb2    = rep("red", 3),
                 xlim      = c(0,1),
                 ylim      = NULL,
                 MAIN      = NULL
                 ){
      ## Load the correct data of volunteer and session
      A <- A[[ V[1], S[1] ]][[ SIDE[1] ]]$data
      B <- B[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
      Aalign <- Aalign[[V[1], S[1]]][[ SIDE[1] ]]
      Balign <- Balign[[V[2], S[2]]][[ SIDE[2] ]]
      ## For the plots
      Ses <- c("A", "B", "C", "D", "E", "F")
      ## Basic informations of the data    
      N1    <- length( A )
      N2    <- length( B )
      maxN  <- max( N1, N2 )
      # basic time grid
      t  <- seq( 0, 1, length.out=100 )
      # arrays for the data
      DATA1 <- DATA2 <- list()
      
      ## Evaluate the geodesic interpolated data either on a constant grid or
      ## on the time registered grid
      if(warping == TRUE){
        for( i in 1:maxN ){
          if( i <= N1 ){
            times <- Aalign$gamma[i,]
            DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Euler")
          }
          if( i <= N2 ){
            times <- Balign$gamma[i,]
            DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Euler")
          }
        }
      }else{
        for( i in 1:maxN ){
          if( i <= N1 ){
            times <- t
            DATA1[[i]] <- eval.geodInterpolation( A[[i]], times = times, out="Euler")  
          }
          if( i <= N2 ){
            times <- t
            DATA2[[i]] <- eval.geodInterpolation( B[[i]], times = times, out="Euler")  
          }
        }
      }
      
      ## Estimate the means of the sessions
      if( MEAN == "Euler" ){
        # initialize means
        mean1 <- mean2 <- c(0,0,0)
        # calculate mean in the chart
        for( i in 1:maxN ){
          if( i <= N1 ){
            mean1 <- mean1 + DATA1[[i]]
          }
          if( i <= N2 ){
            mean2 <- mean2 + DATA2[[i]]
          }
        }
        mean1 <- apply( mean1 / N1, 2, Euler2Quat )
        mean2 <- apply( mean2 / N2, 2, Euler2Quat )
      }else{
        ## Get data in quaternion representation
        DATA1quat <- lapply(DATA1, function(list) apply(list,2, Euler2Quat) )
        DATA2quat <- lapply(DATA2, function(list) apply(list,2, Euler2Quat) )
        ## Compute Ziezold mean on the Sphere
        mean1 <- apply(do.call(rbind, DATA1quat), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
        mean2 <- apply(do.call(rbind, DATA2quat), 2,  function(x) ProjectiveMean( matrix(x,nrow=4))$mean )
      }
      
      ## Get aligning isometry if alignAB == TRUE
      if( alignAB == TRUE ){
        R <- RotEstim( mean1, mean2 )$R
      }else{
        R <- diag(1,4)
      }
      
      ## Get time warping between sessions if warpingAB == TRUE
      if( warpingAB == TRUE ){
        mean1.geod <- geodInterpolation( apply(mean1, 2, Quat2Euler), t  )
        mean2.geod <- geodInterpolation( apply(R %*% mean2, 2, Quat2Euler), t  )
        timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod )
      }
      
      ## Get the new time warped and spatially aligned data2
      for( i in 1:N2 ){
        if(warpingAB == TRUE){
          times <- timeAB$opt.times  
        }else{
          times <- t
        }
        if( warping == TRUE ){
          times <- LinGam( gam = Balign$gamma[i,], grid = t, times = times )
        }
        DATA2[[i]] <- apply( R %*% eval.geodInterpolation( B[[i]], times = times, out="Quat"), 2, Quat2Euler )
      }
      
      ## Confidance bands for volunteer 1
      Angles1 <- GKF.up1 <- GKF.lo1 <- var1 <- list()
      for( ang in 1:3 ){
        Angles1[[ ang ]] <- do.call( cbind, lapply( DATA1, function( list ) list[ ang, ]*Radian ) )
      }
      # mean in Euler chart
      mean1  <- lapply(Angles1, rowMeans )
      # variance in Euler chart
      for( ang in 1 : 3 ){
        var1[[ ang ]]   <- sqrt(rowMeans((Angles1[[ ang ]] - mean1[[ang]])^2) * N1 / (N1-1))
      }
      # GKF construction of simultaneous confidence bands
      for( ang in 1 : 3 ){
        Alpha1  <- gaussianKinematicFormulaT( t( (Angles1[[ ang ]] - mean1[[ ang ]]) / var1[[ang]] ), alpha = alpha / 3 )$c.alpha
        
        GKF.up1[[ ang ]] <- mean1[[ ang ]] + Alpha1*var1[[ ang ]] / sqrt(N1)
        GKF.lo1[[ ang ]] <- mean1[[ ang ]] - Alpha1*var1[[ ang ]] / sqrt(N1)
      }
      # Glue bands of all angles together
      GKF.up1 <- do.call( rbind, GKF.up1 )
      GKF.lo1 <- do.call( rbind, GKF.lo1 )
      mean1   <- do.call( rbind, mean1 )
      
      ## Confidance bands for volunteer 2
      Angles2 <- GKF.up2 <- GKF.lo2 <- var2 <- list()
      for( ang in 1:3 ){
        Angles2[[ ang ]] <- do.call( cbind, lapply( DATA2, function( list ) list[ ang, ]*Radian ) )
      }
      # mean in Euler chart
      mean2  <- lapply(Angles2, rowMeans )
      # variance in Euler chart
      for( ang in 1 : 3 ){
        var2[[ ang ]]   <- sqrt(rowMeans((Angles2[[ ang ]] - mean2[[ang]])^2) * N2 / (N2-1))
      }
      # GKF construction of simultaneous confidence bands
      for( ang in 1 : 3 ){
        Alpha2  <- gaussianKinematicFormulaT( t( (Angles2[[ ang ]] - mean2[[ ang ]]) / var2[[ang]] ), alpha = alpha / 3 )$c.alpha
        
        GKF.up2[[ ang ]] <- mean2[[ ang ]] + Alpha2*var2[[ ang ]] / sqrt(N2)
        GKF.lo2[[ ang ]] <- mean2[[ ang ]] - Alpha2*var2[[ ang ]] / sqrt(N2)
      }
      # Glue bands of all angles together
      GKF.up2 <- do.call( rbind, GKF.up2 )
      GKF.lo2 <- do.call( rbind, GKF.lo2 )
      mean2   <- do.call( rbind, mean2 )
      
      ## Result of the test: If result == 1 there is a significant diference
      ## between the two sessions
      result <- ifelse( any(GKF.lo2 > GKF.up1) || any(GKF.up2 < GKF.lo1), 1, 0 )
      
      ## Show the plot of the test
      if( show.plot == TRUE ){
        ## name  vector for labels in plot
        Angles <- c( "x", "y", "z" )
        ## plot the considered angles
        for( ang in angles){
          if( is.null(ylim) ){
            y.max <- max( c(DATA1[[1]][ang,], DATA2[[1]][ang,]) ) * Radian
            y.min <- min( c(DATA1[[1]][ang,], DATA2[[1]][ang,]) ) * Radian
          }
          if(ang==1){
            MAIN2 <- MAIN
          }else{
            MAIN2 = NULL
          }

          plot(NULL,
               xlim = c( 0, 1 ),
               ylim = c( y.min - 7.5, y.max + 7.5 ),
               main = MAIN2,
               xlab = "standardized frames",
               ylab = paste( Angles[ang],"angle in degree"),
               cex      = 1.5,
               cex.axis = 1.5,
               cex.lab  = 1.5,
               cex.main = 1.5
          )
          lapply(DATA1, function(l) lines( t, l[ang,] * Radian, col=COL1[ang]) )
          lapply(DATA2, function(l) lines( t, l[ang,] * Radian, col=COL2[ang]) )
          if( confbands == TRUE ){
            lines( seq(0,1,length.out=100), GKF.up1[ang,], col=COLcb1, lwd=1 )
            lines( seq(0,1,length.out=100), GKF.lo1[ang,], col=COLcb1, lwd=1 )
            lines( seq(0,1,length.out=100), GKF.up2[ang,], col=COLcb2, lwd=1 )
            lines( seq(0,1,length.out=100), GKF.lo2[ang,], col=COLcb2, lwd=1 )
          }
          if(ang==2){
          legend( "top", paste("Vol", c(V[1],V[2]),
                               "Ses", c(Ses[S[1]],Ses[S[2]]) ),
                  col=c( COL1[ang],COL2[ang] ),
                  lty=rep(1,2),
                  cex=1.4,
                  bty="n" )
          }
        }
      }
 
  ## List with the output    
  list( result = result,
        meanA  = mean1, meanB = mean2,
        cfb.upA = GKF.up1, cfb.loA = GKF.lo1,
        cfb.upB = GKF.up2, cfb.loB = GKF.lo2,
        loBgreaterupA = which(t((GKF.lo2 > GKF.up1)==TRUE)), upBsmallerloB = which(t((GKF.up2 < GKF.lo1)==TRUE)), R=R )
}


TestAB.expModel <- function(
  dataA, dataB, Aalign, Balign,
  V, S, SIDE,
  times     = seq( 0, 1, length.out=100 ),
  angles    = c(1,2,3),
  alpha     = 0.05,
  alignAB   = TRUE,
  warping   = FALSE,
  warpingAB = TRUE,
  show.plot = TRUE,
  colA      = rep("blue", 3), 
  colB      = rep("red", 3),
  xlim      = c(0,1),
  ylim      = NULL,
  MAIN      = NULL
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
    mean2.geod <- geodInterpolation( apply(R %*% mean2, 2, Quat2Euler), times  )
    timeAB     <- timeWarpingSO3( mean1.geod, mean2.geod )
    mean2      <- eval.geodInterpolation( mean2.geod, times = timeAB$opt.times, out="Quat") 
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
  
  #### Test of the same mean
  ## residuals of volunteer 1 with respect to the mean of volunteer 2 and vice versa
  res1 <- Exp.residuals(DATA1, mean2)
  res2 <- Exp.residuals(DATA2, mean1)
  
  ## test whether 0 is included
  result1 <- result2 <- rep(NA, 3)
  up1 <- up2 <- low1 <- low2 <- matrix(NA, nrow=3, ncol=length(times))

  for(a in 1:3){
      ang1     <- do.call(rbind, lapply( res1, function(list) t(list[a,]) ))
      mean1a    <- colMeans(ang1)
      var1     <- sqrt(colMeans(t( (mean1a - t(ang1))^2 ) * N1 / (N1-1)))
      c.alpha1 <- gaussianKinematicFormulaT(ang1 - mean1a, alpha = alpha/6)$c.alpha
      
      ang2     <- do.call(rbind, lapply( res2, function(list) t(list[a,]) ))
      mean2a   <- colMeans(ang2)
      var2     <- sqrt(colMeans(t( (mean2a - t(ang2))^2 ) * N2 / (N2-1)))
      c.alpha2 <- gaussianKinematicFormulaT(ang2 - mean2a, alpha = alpha/6)$c.alpha
      
      up1[a,]  <- mean1a + c.alpha1 * var1 / sqrt(N1)
      low1[a,] <- mean1a - c.alpha1 * var1 / sqrt(N1)
      up2[a,]  <- mean2a + c.alpha2 * var2 / sqrt(N2)
      low2[a,] <- mean2a - c.alpha2 * var2 / sqrt(N2)
      
      if(show.plot==TRUE){
        ## name  vector for labels in plot
        Angles <- c( "x", "y", "z" )
        ## plot the considered angles
          if( is.null(ylim) ){
            y.max <- max( c(up1[a,], up2[a,]) )
            y.min <- min( c(low1[a,], low2[a,]) )
          }
          
          plot(NULL,
               xlim = c( 0, 1 ),
               ylim = c( 1.1*y.min, 1.1*y.max ),
               main = MAIN,
               xlab = "standardized frames",
               ylab = paste( Angles[a],"angle in degree"),
               cex      = 1.5,
               cex.axis = 1.5,
               cex.lab  = 1.5,
               cex.main = 1.5
          )
        lines(times, mean1a, col=colA[a])
        matlines(times, cbind(up1[a,],low1[a,]), col=colA[a],lty=2)
        lines(times, mean2a, col=colB[a])    
        matlines(times, cbind(up2[a,],low2[a,]), col=colB[a],lty=2)
        abline(h=0)
      }
      result1[a] <- !(any(up1[a,]<0) | any(low1[a,]>0))
      result2[a] <- !(any(up2[a,]<0) | any(low2[a,]>0))
  }
  
  ## List with the output    
  list( result1 = result1,
        result2 = result2,
        meanA  = mean1, meanB = mean2,
        upA = up1, lowA = low1,
        upB = up2, lowB = low2
        )
}

