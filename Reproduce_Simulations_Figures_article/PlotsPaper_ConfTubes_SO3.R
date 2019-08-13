#################################### Plots Article #############################
# 
#  This file generates the Plots from the article "..."
#  In order to use it please change the working directory to the directory you
#  saved the git repository.
#
################################################################################

require(KneeMotionAnalytics)
rm(list = ls())     # clear workspace

#################################### Set Constants #############################
WD  <- "/home/drtea/Research/Rpackages/KneeAnally/Reproduce_Simulations_Figures_article" # Change to your folder

setwd(WD)
## Load Wokspace
load("30Apr2015ResultsFullWalk.RData")


###################################### Plots ###################################
speed  = "Walk"
AA      = allGeodWalk
AAalign = GeodWalkSO3align
BB      = allGeodWalk
BBalign = GeodWalkSO3align

times <- seq( 0, 1, length.out=100 )
alpha = 0.1

############### Figure 1a)
V    = c(4,4)
S    = c(3,4)
side = 1

DATA1 <- AA[[V[1],S[1]]][[side]]
DATA2 <- AA[[V[2],S[2]]][[side]]
##  c( bottom, left, top, right)

pdfname <- "Vol4raw.pdf"
pdf(pdfname,  title=pdfname, width=8, height=7)
par(oma = c(0.5, 0.5, 0.5, 0.5),
    mar = c(4.5, 4.5, 0.0, 0.0))
plot(NULL,
     xlim = c(0,100),
     ylim = c(-20, 60),
     main = NULL,
     xlab = "percentage of gait cycle",
     ylab = "Euler angles in [°]",
     cex.lab   = 2.2,
     cex.axis  = 2
)
subPart2 <- 1:length(DATA2$data)
subPart1 <- 1:(length(DATA1$data))
ltyp=c(1,4,5)
for(i in subPart1){for(ang in 1:3){
  l <- length( DATA1$data[[i]]$data[1,] )
  lines( 100*(0:(l-1))/(l-1), DATA1$data[[i]]$data[ang,]*Radian, col="darksalmon", lty=ltyp[ang] )
}}

for(i in subPart2){for(ang in 1:3){
  l <- length( DATA2$data[[i]]$data[1,] )
  lines( 100*(0:(l-1))/(l-1), DATA2$data[[i]]$data[ang,]*Radian, col="cadetblue2", lty=ltyp[ang] )
}}
matlines( 100*seq(0,1,length.out=100), t(DATA1$mean$data)*Radian, col="darkred", lwd=c(2.5,2.5,2.5), lty=1 )
matlines( 100*seq(0,1,length.out=100), t(DATA2$mean$data)*Radian, col="darkblue", lwd=c(2.5,2.5,2.5), lty=1)
legend( "top", inset=.05, c(paste("Vol", V, c("Before","After") )), col=c( "darkred", "darkblue" ), lty=rep(1,2), lwd=c(3.5,3.5), bty="o", box.col="white", bg="white", cex=2.2 )
dev.off()

############### Figure 1b)
pdfname <- "Vol4aligned.pdf"
pdf(pdfname,  title=pdfname, width=8, height=7)
par(oma = c(0.5, 0.5, 0.5, 0.5),
    mar = c(4.5, 4.5, 0.0, 0.0))
      warping   = FALSE
      warpingAB = TRUE
      align     = FALSE
      alignAB   = TRUE
      SIDE      = c(1,1)
      factorN2M = 2
      ## Load the correct data of volunteer and session
      A <- AA[[ V[1], S[1]  ]][[ SIDE[1] ]]$data
      B <- BB[[ V[2], S[2]  ]][[ SIDE[2] ]]$data
      Aalign <- AAalign[[V[1], S[1]]][[ SIDE[1] ]]
      Balign <- BBalign[[V[2], S[2]]][[ SIDE[2] ]]
      ## For the plots
      Ses <- c("A", "B", "A", "B", "E", "F")
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
      
      plot(NULL,
           xlim = c(0,100),
           ylim = c(-20, 60),
           main = NULL,
           xlab = "percentage of gait cycle",
           ylab = "Euler angles in [°]",
           cex.lab   = 2.2,
           cex.axis  = 2
      )
      subPart2 <- 1:length(DATA2)
      subPart1 <- 1:(length(DATA1))
      ltyp=c(1,4,5)
      for(i in subPart1){
        l <- length( DATA1[[i]][1,] )
        matlines( 100*(0:(l-1))/(l-1), t(apply(DATA1[[i]],2,Quat2Euler))*Radian, col="darksalmon", lty=ltyp )
      }
      
      for(i in subPart2){for(ang in 1:3){
        l <- length( DATA2[[i]][1,] )
        matlines( 100*(0:(l-1))/(l-1), t(apply(DATA2[[i]],2,Quat2Euler))*Radian, col="cadetblue2", lty=ltyp )
      }}
      matlines( 100*seq(0,1,length.out=100), t(apply(mean1, 2, Quat2Euler )*Radian), col="darkred", lwd=c(2.5,2.5,2.5), lty=1 )
      matlines( 100*seq(0,1,length.out=100), t(apply(mean2, 2, Quat2Euler )*Radian), col="darkblue", lwd=c(2.5,2.5,2.5), lty=1)
      abline(v=c(0,100), lwd=2, col="black", lty=1)
      abline(v=24, lwd=3, col="black", lty=2)
      abline(v=68, lwd=3, col="black", lty=5)
      abline(v=88, lwd=3, col="black", lty=4)
      legend( "top", inset=.05, c(paste("Vol", V, c("Before","After") )), col=c( "darkred", "darkblue" ), lty=rep(1,2), lwd=c(3.5,3.5), bty="o", box.col="white", bg="white", cex=2.2 )
dev.off()

############### Plot figures for detecting Kneeling effect for volunteers
##### Figure 2,3,4, note that they are put into one pdf.
setwd("/media/sf_Linux/Research/Projects/2017_AOAS_SCBforKnees/Drafts/pics")
V     = c(6,6)
SIDE  = c("Left", "Right")
alpha = 0.05

speed   = "Walk"
AA      = allGeodWalk
AAalign = GeodWalkSO3align
BB      = allGeodWalk
BBalign = GeodWalkSO3align

times <- seq( 0, 1, length.out=100 )

for(side in 1:2){
  pdfname <- paste("Side", SIDE[side],"_alpha",100*alpha,"KneeEffect.pdf",sep="")
  pdf( pdfname,  title = pdfname, width = 13, height = 11 )
  par( oma = c(0.5, 0.5, 0.5, 1.0),
       mar = c(4.5, 5.0, 0.5, 0.5) )
  
  par( mfrow=c(2,2) )
  for(v in 1:8){
    V = c(v,v)
      
      ConfBandsTest.expModelHot(
        dataA = AA, dataB = BB, AAalign, BBalign,
        V = V, S = c(5,3), SIDE = c(side,side),
        times      = seq( 0, 1, length.out=100 ),
        alpha      = alpha,
        show.plot  = TRUE,
        show.plot2 = FALSE,
        Snames     = c("C", "D", "A", "B", "C", "D"),
        xlab       = "Percentage of gait cycle",
        ylab       = "Euler Angles [°]",
        cex.lab    = 2.5,
        cex.axis   = 1.5
      )
      
      ConfBandsTest.expModelHot(
        dataA=AA, dataB=BB, AAalign, BBalign,
        V=V, S=c(5,4), SIDE=c(side,side),
        times     = seq( 0, 1, length.out=100 ),
        alpha     = alpha,
        show.plot = TRUE,
        show.plot2 = FALSE,
        Snames    = c("C", "D", "A", "B", "C", "D"),
        xlab       = "Percentage of gait cycle",
        ylab       = "Euler Angles [°]",
        cex.lab    = 2.5,
        cex.axis   = 1.5
      )
      
      ConfBandsTest.expModelHot(
        dataA=AA, dataB=BB, AAalign, BBalign,
        V=V, S=c(6,3), SIDE=c(side,side),
        times     = seq( 0, 1, length.out=100 ),
        alpha     = alpha,
        show.plot = TRUE,
        show.plot2 = FALSE,
        Snames    = c("C", "D", "A", "B", "C", "D"),
        xlab       = "Percentage of gait cycle",
        ylab       = "Euler Angles [°]",
        cex.lab    = 2.5,
        cex.axis   = 1.5
      )
      
      
      ConfBandsTest.expModelHot(
        dataA = AA, dataB = BB, AAalign, BBalign,
        V = V, S = c(6,4), SIDE = c(side,side),
        times      = seq( 0, 1, length.out=100 ),
        alpha      = alpha,
        show.plot  = TRUE,
        show.plot2 = FALSE,
        Snames     = c("C", "D", "A", "B", "C", "D"),
        xlab       = "Percentage of gait cycle",
        ylab       = "Euler Angles [°]",
        cex.lab    = 2.5,
        cex.axis   = 1.5
      )
  }
  dev.off()  
  }

##### Plot of SCB test for non kneeling sessions, note that they are put into one pdf.
for(side in 1:2){
  pdfname <- paste("Vols_Side", SIDE[side],"_alpha",100*alpha,".pdf",sep="")
  pdf( pdfname,  title = pdfname, width = 13, height = 11 )
  par( oma = c(0.5, 0.5, 0.5, 1.0),
       mar = c(4.5, 5.0, 0.5, 0.5) )
  
    for(v in 1:8){
    V = c(v,v)
    
    ConfBandsTest.expModelHot(
      dataA = AA, dataB = BB, AAalign, BBalign,
      V = V, S = c(3,4), SIDE = c(side,side),
      times      = seq( 0, 1, length.out=100 ),
      alpha      = alpha,
      show.plot  = TRUE,
      show.plot2 = FALSE,
      Snames     = c("C", "D", "A", "B", "C", "D"),
      xlab       = "Percentage of gait cycle",
      ylab       = "Euler Angles [°]",
      cex.lab    = 2.5,
      cex.axis   = 1.5
    )
    }
  dev.off()
}

############################ Generate table 2 ############################
Trials = matrix(NA, 8, 4)

for( v in 1:8 ){for( s in 3:6 ){
  Trials[v,s-2] = length(allGeodWalk[[v,s]][[1]]$data)
}}

MedTrials  = apply(Trials, 2, quantile)
colnames(MedTrials) = c("A", "B", "C", "D")
print(xtable(t(MedTrials)))