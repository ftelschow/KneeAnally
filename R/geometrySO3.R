################################################################################
#
#     Functions computing the geometry of SO(3)
#
################################################################################
#' Computes the Lie Exponential of SO(3)
#'
#' @param v a vector
#' @return numeric p-norm of the vector
#' @export
#'
vecNorm <- function( v, p=2 ){
  ( sum(v^p))^(1/p)
}

#' Computes the Lie Exponential of SO(3)
#'
#' @param A skew symmetric matrix or vector
#' @return Matrix rotation matrix
#' @export
#'
Exp.SO3 <- function(A){
  ## Transform vector input into scew matrix
  if( is.vector(A) ){
    A <- Vec2screwM(A)
  }
  ## Chirikjian (2000) p.120
  if( sqrt( sum(diag(t(A)%*%A))/2 ) != 0 ){
        absA <- sqrt(sum(diag(t(A)%*%A)) / 2)
        R <- diag(rep(1,3)) + sin(absA) / absA * A + (1-cos(absA)) / absA^2 * A%*%A
  }else{
        R <- diag(rep(1,3))
  }
  R
}

#' Computes the log map, i.e. the inverse of the Lie Exponential of SO(3)
#'
#' @param R rotation matrix
#' @return Matrix skew symmmetric matrix
#' @export
#'
Log.SO3 <- function(R, eps=1e-8){
  t = sum(diag(R))
  if( t >= 3-eps ){
    ang = acos( (t-1)/2 )
    X <- ( 0.5 + ang^2/12 ) * (R-t(R))
  }else if (t <= -1+eps ){
    q <- Rot2Quat(R)
    ang = sqrt( sum(q[2:4]^2) )
    X <- 2*asin(ang) * Vec2screwM(q[2:4]) / ang
  }else{
    theta <- acos( (t-1)/2 )
    X <- (theta * R - theta * t(R)) / ( 2*sin(theta) )
  }
  X
}

#' Computes the intrinsical distance between two elements of S^n w.r.t. the
#' Fubini Study metric
#'
#' @param q a unit vector
#' @param p a unit vector
#' @return Matrix skew symmmetric matrix
#' @export
#'
IntrinsicDist.SO3 <- function(q,p){
  if( !is.vector(q) ){
    q = Rot2Quat(q)
  }
  if( !is.vector(p) ){
    p = Rot2Quat(p)    
  }
  
  if(sum(q*p)>1){
    phi = 0
  } else if( sum(q*p)< -1){
    phi = pi
  }else{
    ## angle between q and p (<pi)
    phi <- acos(sum(q*p))
  }
  ## return length between q and p
  2*min(phi,pi-phi)
}

#' Computes the standard inner product between two elements of the Lie algebra
#' so(3) or more generally on all matrices
#'
#' @param X a 
#' @param Y a unit vector
#' @return Number value of the inner product
#' @export
#'
ska.SO3 <- function(X,Y){
  1/2*sum(diag(t(Y)%*%X))
}

#' Computes the angle between two matrices
#'
#' @param X a 
#' @param Y a unit vector
#' @return Number angle between the matrices
#' @export
#'
angle.SO3 <- function(X,Y){
  acos( ska.SO3(X,Y) / sqrt(ska.SO3(X,X)*ska.SO3(Y,Y)) )
}

#' Maps vector two skew symmetric matrix
#'
#' @param x a vector
#' @return Matrix skew symmetric matrix
#' @export
#'
Vec2screwM <- function(x){
  
  cbind(c(0,x[3], -x[2]), c(-x[3],0,x[1]), c(x[2],-x[1],0))
  
}
#' Maps screw matrix into vector representation
#'
#' @param X a Matrix
#' @return Vector
#' @export
#'
screwM2Vec <- function(X){
  c( X[3,2] , X[1,3] ,X[2,1] )
}

#' Computes the total variation distance of curves in SO(3) 
#'
#' @param f list geodesicInterpolation object of a curve.
#' @param g list geodesicInterpolation object of a curve.
#' @return numeric TV distance between f and g
#' @export
#'
wFun <- function( f, k, g, n, m ) {
  w1 = g[,n] ; w2 = g[,m] ;
  v1 = f[,k] ; v2 = f[,k+1] ;
  
  p1 = QuatMult(w1, QuatInv(v1))
  p2 = QuatMult(w2, QuatInv(v2))
  p1 = p1 / sqrt( sum(p1*p1) )
  p2 = p2 / sqrt( sum(p2*p2) )
  
  IntrinsicDist.SO3( p1, p2 )
}

#' Computes the total variation distance of curves in SO(3) 
#'
#' @param f list geodesicInterpolation object of a curve.
#' @param g list geodesicInterpolation object of a curve.
#' @return numeric TV distance between f and g
#' @export
#'
length.curve <- function(f, g){
  N <- dim(f)[2]
  
  sum(vapply(1:(N-1), function(n) wFun( f, n, g, n, n+1), FUN.VALUE=0 ) )
}

#' Computes the total variation distance of curves in SO(3) 
#'
#' @param f list geodesicInterpolation object of a curve.
#' @param g list geodesicInterpolation object of a curve.
#' @return numeric TV distance between f and g
#' @export
#'
dist.FT1 <- function( f, g ) {
  N <- dim(f)[2]
  p = vapply(1:N, function(n) QuatMult( g[,n], QuatInv(f[,n])), rep(1,4)  )
  p = apply(p, 2, function(n) n/sqrt(sum(n^2)))

  sum( vapply(1:(N-1), function(n) IntrinsicDist.SO3( p[,n], p[,n+1]), FUN.VALUE=1) )
}
#' @export
#'
dist.FT2 <- function( f, g ) {
  N <- dim(f)[2]
  p = vapply(1:N, function(n) QuatMult( QuatInv(g[,n]), f[,n]), rep(1,4)  )
  p = apply(p, 2, function(n) n/sqrt(sum(n^2)))
  
  sum( vapply(1:(N-1), function(n) IntrinsicDist.SO3( p[,n], p[,n+1]), FUN.VALUE=1) )
}
#' @export
#'
dist.FT <- function( f, g ) {
  (dist.FT1(f,g) + dist.FT2(f,g)) / 2
}

#' @export
dist.L2 <- function( f, g ) {
  N <- dim(f)[2]
  sum(vapply(1:N, function(n) IntrinsicDist.SO3( f[,n], g[,n])^2, FUN.VALUE=1  ) )
}

#' @export
minDistL2 <- function(p, m1, m2){
  P   <- Exp.SO3( Vec2screwM(p[1:3]) )
  m2p <- apply(m2, 2, function(col) Rot2Quat(P%*%Quat2RotC(col)%*%t(P) ) )
  
  dist.L2( m,m2p )
}

#' @export
hatPillQ <- function(p, m1, m2){
  Q   <- ExpSO3C( Vec2screwMC(p[1],p[2],p[3]) )
  m2p <- apply(m2, 2, function(col) Rot2QuatC(Quat2RotC(col)%*%Q ) )
  
  distFT1C( m1,m2p )
}

#' @export
hatPillP <- function(p, m1, m2){
  P   <- ExpSO3C( Vec2screwMC(p[1],p[2],p[3]) )
  m2p <- apply(m2, 2, function(col) Rot2QuatC(P%*%Quat2RotC(col) ) )
  
  distFT1C( apply(m1,2,QuatInvC), apply(m2p,2,QuatInvC) )
}