################################################################################
#     Functions for geodesic interpolation on SO(3)
#
################################################################################
#' Geodesic Interpolation on SO(3) 
#'
#' @param DATA Matrix of Euler angles 3 x N or quaternions 4 x N
#' @param times Vector of measurement times
#' @return list with elements
#'  \itemize{
#'   \item data DATA Matrix of Euler angles 3 x N
#'   \item directions list of symmetric matrices, i.e. geodesic directions
#'   \item times Vector of measurement times
#' }
#' @export
#' 
geodInterpolation <- function(DATA, times){
  ## length of data
  N <- length(times)
  
  ## direction vector of the geodesic
  geod.direction <- list()
    
    for( n in 1:(N-1) ){
      ## Get geodesic direction
      if( length(DATA[,1]) == 3 ){
          rotk                    <- Euler2Rot(DATA[,n])
          rotk1                   <- Euler2Rot(DATA[,n+1])
      }else if( length(DATA[,1]) == 4 ){
        rotk                    <- Quat2Rot(DATA[,n])
        rotk1                   <- Quat2Rot(DATA[,n+1])        
      }else{
        print("Wrong Input data: Only 3xN or 4xN matrices allowed")
      }
      geod.direction[[n]]     <- LogSO3C(t(rotk)%*%rotk1)
    }
  
    if( length(DATA[,1]) == 4 ){
      DATA <- apply(DATA, 2, Quat2Euler )
    }
  
  list( data = DATA, directions = geod.direction, times = times )  
}

#' Evaluate Geodesic Interpolation on SO(3) 
#'
#' @param InterPol List created by geodInterpolation() 
#' @param times Vector where to evaluate the interpolation
#' @param out string specifying output: either "Euler" or "Quat"
#' @return gamma Matrix (3 x N) of SO(3) interpolated euler angles
#' @export
#' 
eval.geodInterpolation <- function(InterPol, times, out="Euler"){
  ## time grid of the interpolation
  measurement.times <- InterPol$times
  ## Output vector
  if(out=="Euler"){
    gamma.new <- matrix(NA,3,length(times))
  }else{
    gamma.new <- matrix(NA,4,length(times))   
  }

  ## Compute interpolated values
  for(n in 1:length(times)){
    ## On which geodesic segment does it lie?
    if(times[n]==measurement.times[1]){
      k=2
    }else{
      k <- which(times[n] <= measurement.times)[1]  
    }
    ## fraction on the geodesic segment
    dt <- (times[n]-measurement.times[k-1]) /
                                  (measurement.times[k]-measurement.times[k-1])
    ## Compute SO(3) value of new time points
    if(out=="Euler"){
            gamma.new[,n] <- Rot2Euler(
                         Euler2Rot(InterPol$data[,k-1])%*%
                                   ExpSO3C(dt*InterPol$dir[[k-1]])
                      )
    }else{
      gamma.new[,n] <- Rot2QuatC(Euler2Rot(InterPol$data[,k-1])%*%
                                  ExpSO3C(dt*InterPol$dir[[k-1]]))      
    }
  }
  gamma.new
}

#' Inverse of a piecwise linear function f.
#'
#' @param y Numeric evaluation point of inverse.
#' @param fx Vector containing the values of f at the base points.
#' @param times Vector base times corresponding to the values of fx.
#' @return val Numeric value of the inverse at y.
#' @export
#' 
LinInterpolationInverse <- function(y, fx, times){
  val <- 0
  if(y < max(times)){
  n <- which( times > y )[1] - 1
  val <-  diff(times)[n] / diff(fx)[n] * ( y - fx[n]) + times[n]
  } else {
    n <- length(times) - 1
    val <- diff(times)[n] / diff(fx)[n] * ( y - fx[n]) + times[n]    
  }
  val
}