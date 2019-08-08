################################################################################
#
#     Representations of SO(3)
#
################################################################################
#' Rotation matrix about a coordinate axis
#'
#' @param angle Number rotation angle in radian
#' @param coord Char either "x","y" or "z", specifying the coordinate axis of
#' the rotation
#' @return Matrix Rotation matrix around the specidied coordinate axis
#' @export
#' 
Euler2RotMatrix <- function( angle, coord="x" ) {
  R <- NULL
  if (coord == "x"){
    R <- cbind(c(1,0,0),c(0,cos(angle),sin(angle)),c(0,-sin(angle),cos(angle)))
  }else if (coord == "y"){ 
    R <- cbind(c(cos(angle),0,-sin(angle)),c(0,1,0),c(sin(angle),0,cos(angle)))
  }else if (coord == "z"){ 
    R <- cbind(c(cos(angle),sin(angle),0),c(-sin(angle),cos(angle),0),c(0,0,1))
  }else{
    print("Error wrong Input for coord: Please choose either x,y or z")
  }
  R
}

#' Quaternion q representing a rotation matrix about a coordinate axis via
#' Rv = Adj(q,c(0,v))
#'
#' @param angle Number rotation angle in radian
#' @param coord Char either "x","y" or "z", specifying the coordinate axis of
#' the rotation
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
Euler2RotQuat <- function(angle,coord="x") {
  if (coord == "x") 
    q <- c(cos(angle/2), c(sin(angle/2),0,0))
  if (coord == "y") 
    q <- c(cos(angle/2), c(0,sin(angle/2),0))
  if (coord == "z") 
    q <- c(cos(angle/2), c(0,0,sin(angle/2)))
  q
}

#' Euler angles to quaternion
#'
#' @param angle Vector three euler angles in radian
#' @param convention Choose an Euler angle convention currently "yxz" and "zxy"
#' @return Vector with for elements, i.e. quaternion q representing the rotation
#' @export
#' 
Euler2Quat <- function(angles, convention="yxz") {
  # if( convention = "yxz" ){
  #   q = QuatMultC( Euler2RotQuat(angles[2],"y"),QuatMultC(Euler2RotQuat(angles[1],"x"), Euler2RotQuat(angles[3],"z")) )
  # }else if( convention = "zxy" ){
  #   q = QuatMultC( Euler2RotQuat(angles[3],"z"),QuatMultC(Euler2RotQuat(angles[1],"x"), Euler2RotQuat(angles[2],"y")) )
  # }
  # q
  Rot2Quat( Euler2Rot(angles) )
}

#' Euler angles to rotation matrix
#'
#' @param angle Vector three euler angles in radian
#' @param convention Choose an Euler angle convention currently "yxz" and "zxy"
#' @return Matrix 3x3 rotation matrix
#' @export
#' 
Euler2Rot <- function( angles, convention = "yxz" ) {
  if( convention == "yxz" ){
    R = Euler2RotMatrix(angles[2],"y")%*%Euler2RotMatrix(angles[1],"x")%*% Euler2RotMatrix(angles[3],"z")
  }else if( convention == "zxy" ){
    R = Euler2RotMatrix(angles[3],"z")%*%Euler2RotMatrix(angles[1],"x")%*% Euler2RotMatrix(angles[2],"y")
  }
  R
}

#' Quaternion to rotation matrix
#'
#' @param q vector with 4 elements, i.e. a quaternion
#' @return Matrix 3x3 rotation matrix
#' @export
#' 
Quat2Rot <- function(q) {
  a <- q[1]; b <- q[2]; c <- q[3]; d <- q[4];
  cbind(
        c(  1-2*c^2-2*d^2, 2*(b*c+a*d)     , 2*(b*d-a*c)   ),
        c(  2*(b*c-a*d)  , 1-2*b^2-2*d^2   , 2*(c*d+a*b)   ),  
        c(  2*(b*d+a*c)  , 2*(c*d-a*b)     , 1-2*b^2-2*c^2 )
      )
}

#' Quaternion to Euler angles
#'
#' @param q vector with 4 elements, i.e. a quaternion
#' @return vector three Euler angles
#' @export
#' 
Quat2Euler <- function(q) {
  Rot2Euler(Quat2Rot(q))
}

#' Matrix to Euler angles
#'
#' @param R rotation matrix
#' @return vector three Euler angles
#' @export
#' 
Rot2Euler <- function(R) {
  if( abs(R[2,3]) < 1 ){
    x <- asin(-R[2,3])
    z <- atan2( R[2,1] / cos(x) , R[2,2] / cos(x) )
    y <- atan2( R[1,3] / cos(x) , R[3,3] / cos(x) )
  }else{
    
  }
  euler <- rep(NA,3)
  names(euler) <- c("x","y","z")
  euler <- c(x,y,z)
  euler
}

#' #' Matrix to quaternion
#' #'
#' #' @param R
#' #' @return vector with 4 elements, i.e. a quaternion
#' #' @export
#' #' 
#' Rot2Quat <- function(R) {
#'   s <- sqrt( sum(diag(R))+1 ) / 2
#'   if( s !=0 ){
#'     q =  c( s, c( R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2] ) / (4*s) )
#'   }else if( R[1,2]==0 & R[2,3] == 0 ){
#'     q = c( sqrt( (1+R[1,1])/2 ) , 0, sqrt( (1-R[1,1])/2 ),0 )
#'   }else if( R[1,3]==0 & R[2,3] == 0 ){
#'     q = c( sqrt( (1+R[1,1])/2 ) , 0, 0, sqrt( (1-R[1,1])/2 ) )
#'   }else if( R[1,2]==0 & R[2,3] == 0 ){
#'     q =  c( sqrt( (1+R[1,1])/2 ), sqrt( (1-R[1,1])/2 ), 0, 0 )
#'   }else{
#'     q = c( 0, R[1,2]*R[1,3], R[1,2]*R[2,3], R[1,3]*R[2,3] )
#'   }
#'   q / sqrt(sum(q^2))
#' }

#' #' Matrix to quaternion
#'
#' @param R
#' @return vector with 4 elements, i.e. a quaternion
#' @export
#'
Rot2Quat <- function(R) {
  if( sum(diag(R))==3 ){
    q = c(1, 0, 0, 0)
  }else{
    # find biggest diagonal entry for stable computation
    divisor = c( sqrt(abs(1+sum(diag(R))))/2, sqrt(abs(1+sum(diag(R)*c(1,-1,-1))))/2, sqrt(abs(1+sum(diag(R)*c(-1,1,-1))))/2, sqrt(abs(1+sum(diag(R)*c(-1,-1,1))))/2 )
    m = which.max(divisor)
    q = rep(NA,4)
    
    if( m == 1 ){
      q[1] <- divisor[m]
      q[2] <- ( R[3,2] - R[2,3] ) / ( 4*divisor[m] )
      q[3] <- ( R[1,3] - R[3,1] ) / ( 4*divisor[m] )
      q[4] <- ( R[2,1] - R[1,2] ) / ( 4*divisor[m] )
    }else if( m == 2){
      q[1] <- ( R[3,2] - R[2,3] ) / ( 4*divisor[m] )
      q[2] <- divisor[m]
      q[3] <- ( R[1,2] + R[2,1] ) / ( 4*divisor[m] )
      q[4] <- ( R[1,3] + R[3,1] ) / ( 4*divisor[m] )
    }else if( m == 3 ){
      q[1] <- ( R[1,3] - R[3,1] ) / ( 4*divisor[m] )
      q[2] <- ( R[1,2] + R[2,1] ) / ( 4*divisor[m] )
      q[3] <- divisor[m]
      q[4] <- ( R[2,3] + R[3,2] ) / ( 4*divisor[m] )
    }else{
      q[1] <- ( R[2,1] - R[1,2] ) / ( 4*divisor[m] )
      q[2] <- ( R[1,3] + R[3,1] ) / ( 4*divisor[m] )
      q[3] <- ( R[2,3] + R[3,2] ) / ( 4*divisor[m] )
      q[4] <- divisor[m]
    }
  }
  # map q to the quaternion with positiv first entry
  if( q[1]<0 ){
    q=-q
  }
  q / sqrt(sum(q^2))
}

#' SO(4) Matrix to two quaternions
#'
#' @param R
#' @return matrix  4x2 elements, containing q_l and q_R
#' @export
#'
Rot42Quat <- function(R) {
  # Compute associate matrix
  M = 1/4 * rbind(
  c( R[1,1]+R[2,2]+R[3,3]+R[4,4], R[2,1]-R[1,2]-R[4,3]+R[3,4], R[3,1]+R[4,2]-R[1,3]-R[2,4], R[4,1]-R[3,2]+R[2,3]-R[1,4] ) ,
  c( R[2,1]-R[1,2]+R[4,3]-R[3,4], -R[1,1]-R[2,2]+R[3,3]+R[4,4], R[4,1]-R[3,2]-R[2,3]+R[1,4], -R[3,1]-R[4,2]-R[1,3]-R[2,4] ) ,
  c( R[3,1]-R[4,2]-R[1,3]+R[2,4], -R[4,1]-R[3,2]-R[2,3]-R[1,4], -R[1,1]+R[2,2]-R[3,3]+R[4,4], R[2,1]+R[1,2]-R[4,3]-R[3,4] ) ,
  c( R[4,1]+R[3,2]-R[2,3]-R[1,4], R[3,1]-R[4,2]+R[1,3]-R[2,4], -R[2,1]-R[1,2]-R[4,3]-R[3,4], -R[1,1]+R[2,2]+R[3,3]-R[4,4] )
  )
  SVD=svd(M)
  return( list( p=SVD$u[,1], q=SVD$v[,1] ) )
  
}