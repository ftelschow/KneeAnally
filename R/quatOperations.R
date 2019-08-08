################################################################################
#
#     Operations on quaternions
#
################################################################################
#' Multiplication of two quaternions
#'
#' @param q1 vector with 4 elements, i.e. a quaternion
#' @param q2 vector with 4 elements, i.e. a quaternion
#' @return vector Multiplication of q1 and q2
#' @export
#' 
QuatMult <- function(q1,q2) {
  c( q1[1] * q2[1] - crossprod(q1[2:4], q2[2:4]), q1[1]*q2[2:4]+q2[1]*q1[2:4]+
      q1[c(3,4,2)]*q2[c(4,2,3)]-q2[c(3,4,2)]*q1[c(4,2,3)])
}

#' Conjugation of a quaternion
#'
#' @param q vector with 4 elements, i.e. a quaternion
#' @return vector Conjugation of q
#' @export
#' 
QuatConj <- function(q){ c(q[1],-q[2:4]) }

#' Conjugation of a quaternion
#'
#' @param q vector with 4 elements, i.e. a quaternion
#' @return vector Inverse of q
#' @export
QuatInv <- function(q){ c(q[1],-q[2:4]) / sum(q^2) }

#' Adjungation of a quaternion by another quaternion
#'
#' @param p vector with 4 elements, i.e. a quaternion
#' @param q vector with 4 elements, i.e. a quaternion
#' @return vector Adjungation of q by p. i.e. pqp^{-1}
#' @export
#' 
Adj <- function(p,q) QuatMult(p,QuatMult(q,QuatInv(p)))