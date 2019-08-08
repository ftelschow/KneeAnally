################################################################################
#
#   Ball's functional alignment
#
################################################################################
KreuzProd <- function(a,b){
  c( a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1] )
}

Rot2AxisAngle <-  function(R){
#   v = c( R[3,2] - R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])
#   norm <- sqrt(sum(v^2))
#   theta <- atan(norm / (R[1,1]+R[2,2]+R[3,3]-1))
  v <- screwM2Vec(Log.SO3(R))
  theta <- sqrt(sum(v^2))
  list(u = v / theta, ang = theta)
#    list(u = v / norm, ang = theta)
}

AxisAngle2Rot <- function(AxisAngle){
  u <- AxisAngle$u
  ang <- AxisAngle$ang
  rbind(
        c( u[1]*u[1]*(1-cos(ang)) + cos(ang),      u[1]*u[2]*(1-cos(ang)) - u[3]*sin(ang), u[1]*u[3]*(1-cos(ang)) + u[2]*sin(ang) ),
        c( u[2]*u[1]*(1-cos(ang)) + u[3]*sin(ang), u[2]*u[2]*(1-cos(ang)) + cos(ang), u[2]*u[3]*(1-cos(ang)) - u[1]*sin(ang) ),
        c( u[3]*u[1]*(1-cos(ang)) - u[2]*sin(ang), u[3]*u[2]*(1-cos(ang)) + u[1]*sin(ang), u[3]*u[3]*(1-cos(ang)) + cos(ang) ) 
    )
}

spat.ave <- function(A){
  S = Quat2Rot(A[,1])
  for(i in 2:dim(A)[2]){
    G <- t(S) %*% Quat2Rot(A[,i])
    u <- Rot2AxisAngle(G)
    u$ang <- u$ang / i
    P <- AxisAngle2Rot( u )
    S <- S %*% P
  }
  S
}

BallAlignment <- function(A, r=c(0,1,0), ave = "Ball"){
  ## Compute spatial average
  if( ave == "Ball" ){
    A.ave <- spat.ave(A)
  }else{
    A.ave <- Quat2Rot(ProjectiveMean(A)$mean)
  }
  
  ind.pos <- NULL
  for(i in 1:dim(A)[2]){
    R <- t(A.ave) %*% Quat2Rot(A[,i])
    if( sum(Rot2AxisAngle(R)$u * r) > 0 ){
      ind.pos <- c(ind.pos,i)
    }
  }
  
  if( ave == "Ball" ){
    A.p <- spat.ave( A[,ind.pos] )
    A.n <- spat.ave( A[,setdiff(1:dim(A)[2], ind.pos)] )
  }else{
    A.p <- Quat2Rot(ProjectiveMean( A[,ind.pos] )$mean)
    A.n <- Quat2Rot(ProjectiveMean( A[,setdiff(1:dim(A)[2], ind.pos)] )$mean)
  }
  
  ## Get M.A
  J.M <- t(A.n) %*% A.p
  m   <- Rot2AxisAngle(J.M)$u
  theta.m <- acos( sum(m * r) )
  if(theta.m > pi/2){
    theta.m = theta.m - pi
    m <- -m
  }
  u.m   <- KreuzProd(m,r)
  u.m   <- u.m / sqrt(sum(u.m^2))
  M.A <- AxisAngle2Rot( list(u=u.m, ang=theta.m) )
  
  ## Get B.A
  J.B <- A.n %*% t(A.p)
  b   <- Rot2AxisAngle(J.B)$u
  theta.b <- acos( sum(b * r) )
  if(theta.b > pi/2){
    theta.b = theta.b - pi
    b <- -b
  }
  u.b   <- KreuzProd(b,r)
  u.b   <- u.b / sqrt(sum(u.b^2))
  B.A <- AxisAngle2Rot( list(u=u.b, ang=theta.b) )
  
  list(rR = M.A, lR = t(B.A))
}