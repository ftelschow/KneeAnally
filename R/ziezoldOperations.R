#' Bring elemtents in optimal position w.r.t. group action of S^0 on R^k
#'
#' @param A matrix (k x N) containing N vectors in R^k
#' @return A matrix (k x N) containing N vectors in R^k in optimal position
#' @export
#'
OptPos.Quat <- function(A, j=1){
  if(length(A[1,])<2) return("Wrong Input")
    q.ref <- A[,j]; n <- length(A[1,]);
  if(n>2){
    for(i in 1:n){
      if(i!=j && sum((q.ref-A[,i])^2) > sum((q.ref+A[,i])^2)) A[,i] <- -A[,i]
    }
  }
  A
}
#' Ziezold distance on R^k / S^0
#'
#' @param x vector in R^k
#' @param y vector in R^k
#' @return number value of the distance between x and y
#' @export
#'
Ziez.dist <- function(x,y){
  min( sqrt(sum((x-y)^2)),sqrt(sum((x+y)^2)) )
}


## SO(3)
#' Calculate Ziezold mean on R^k / S^0 for a sample
#'
#' @param L matrix (k x N) with columns containing the realisations of the r.v.
#' @param MaxIt number giving the maximal number of iterations
#' @param err threshold for the iteration
#' @return list with elements
#'  \itemize{
#'   \item x vector containing the ziezold mean of the sample
#'   \item OptPos matrix k x N contains the data in optimal position to x
#'   \item Variances vector containing the variances for each observation
#' }
#' @export
#'
ProjectiveMean <- function(L, MaxIt=100, err=10^(-10)){
  x <- L[,1]
  cont <- TRUE; n <- 0
  while (cont==TRUE) {
    n <- n+1
    x.old  <- x
    OptPos  <- OptPosQuatC(L, x)
    
    dim(OptPos) <- dim(L)
    y      <- apply(OptPos,1,mean); x <- y/sqrt(sum(y^2))
    d1     <- sqrt(sum((x-x.old)^2)); d2<- sqrt(sum((x+x.old)^2))
    if (d1 < err || d2 < err || n > MaxIt) cont <- FALSE
  }
  if (cont==FALSE)
  { 
    list(
          mean=x, optPos=OptPos, variances=apply(
          OptPos,2,function(y) sum((y-x)^2))
         )
  }
  else return("FAIL")
}