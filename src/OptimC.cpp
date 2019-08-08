#include <RcppArmadillo.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;
////////////////////////////////////////////////////////////////
//            Programming tools
////////////////////////////////////////////////////////////////
// NumericMatrix to arma::matrix
// [[Rcpp::export]]
arma::mat NumMat2armaMat (NumericMatrix x){
  arma::mat y = as<arma::mat>(x) ;
  return(y) ;
}
// arma::matrix to NumericMatrix
// [[Rcpp::export]]
NumericMatrix armaMat2NumMat (arma::mat x){
  NumericMatrix y = wrap(x) ;
  return(y) ;
}

////////////////////////////////////////////////////////////////
//            Operations on Quaternions
////////////////////////////////////////////////////////////////
//Function computing quaternion multiplication
// [[Rcpp::export]]
NumericVector QuatMultC ( NumericVector q1, NumericVector q2 ) {
  int n = q1.size() ;
  NumericVector out(n) ;
  
  out(0) = q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2)-q1(3)*q2(3) ;
  out(1) = q1(0)*q2(1) + q1(1)*q2(0) + q1(2)*q2(3)-q1(3)*q2(2) ;
  out(2) = q1(0)*q2(2) + q1(2)*q2(0) + q1(3)*q2(1)-q1(1)*q2(3) ;
  out(3) = q1(0)*q2(3) + q1(3)*q2(0) + q1(1)*q2(2)-q1(2)*q2(1) ;
  
  return( out ) ;
}
//Function computing quaternion multiplication
// [[Rcpp::export]]
arma::vec QuatMultC2 ( arma::vec q1, arma::vec q2 ) {
  int n = q1.size() ;
  arma::vec out(n) ;
  
  out(0) = q1(0)*q2(0) - q1(1)*q2(1) - q1(2)*q2(2)-q1(3)*q2(3) ;
  out(1) = q1(0)*q2(1) + q1(1)*q2(0) + q1(2)*q2(3)-q1(3)*q2(2) ;
  out(2) = q1(0)*q2(2) + q1(2)*q2(0) + q1(3)*q2(1)-q1(1)*q2(3) ;
  out(3) = q1(0)*q2(3) + q1(3)*q2(0) + q1(1)*q2(2)-q1(2)*q2(1) ;
  
  return( out ) ;
}
//Function computing quaternion inverse
// [[Rcpp::export]]
NumericVector QuatInvC ( NumericVector q ) {
  NumericVector out(4) ;
  
  double normq = sum(q*q) ;
  
  out(0) =  q(0) / normq ;
  out(1) = -q(1) / normq ;
  out(2) = -q(2) / normq ;
  out(3) = -q(3) / normq ;
  
  return( out ) ;
}
//Function computing quaternion inverse
// [[Rcpp::export]]
arma::vec QuatInvC2 ( arma::vec q ) {
  arma::vec out(4) ;
  
  double normq = sum(q%q) ;
  
  out(0) =  q(0) / normq ;
  out(1) = -q(1) / normq ;
  out(2) = -q(2) / normq ;
  out(3) = -q(3) / normq ;
  
  return( out ) ;
}

////////////////////////////////////////////////////////////////
//            Conventions on Rotations
////////////////////////////////////////////////////////////////
//Function computing Euler angles from rotation matrix
// [[Rcpp::export]]
NumericVector Rot2EulerC ( arma::mat R) {
  NumericVector e(3);
  
  if( R(1,2) < 1 && R(1,2)>-1 ){
    e(0) = asin(-R(1,2));
    e(2) = atan2( R(1,0) / cos(e(0)) , R(1,1) / cos(e(0)) );
    e(1) = atan2( R(0,2) / cos(e(0)) , R(2,2) / cos(e(0)) );
  }else{
  }
  return( e );
}
// Function computing quaternion from euler angles
// [[Rcpp::export]]
arma::vec Euler2QuatCe ( double angle, CharacterVector coord) {
  arma::vec q(4);
  
  if (LogicalVector(coord(0) == "x")(0)) {
    q(0) = cos(angle/2);
    q(1) = sin(angle/2);
    q(2) = 0;
    q(3) = 0;
  } else if (LogicalVector(coord(0) == "y")(0)) {
    q(0) = cos(angle/2);
    q(1) = 0;
    q(2) = sin(angle/2);
    q(3) = 0;
  } else {
    q(0) = cos(angle/2);
    q(1) = 0;
    q(2) = 0;
    q(3) = sin(angle/2);
  }
  
  return( q );
}
// Function computing quaternion from euler angles
// [[Rcpp::export]]
NumericVector Euler2QuatC ( arma::vec angles ) {
  return( QuatMultC( as<NumericVector>(wrap(Euler2QuatCe( angles(1), CharacterVector::create("y")))), QuatMultC( as<NumericVector>(wrap(Euler2QuatCe( angles(0), CharacterVector::create("x")))),
                                       as<NumericVector>(wrap(Euler2QuatCe( angles(2), CharacterVector::create("z"))))) ) ) ;
}
// Function transforming a quaternion into a rotation matrix
// [[Rcpp::export]]
arma::mat Quat2RotC (arma::vec q) {
  arma::mat R(3,3);
  double a = q(0);
  double b = q(1);
  double c = q(2);
  double d = q(3);
  
  R(0,0) = 1-2*c*c-2*d*d; R(1,0) = 2*(b*c+a*d)     ; R(2,0) = 2*(b*d-a*c)   ;
  R(0,1) = 2*(b*c-a*d)  ; R(1,1) = 1-2*b*b-2*d*d   ; R(2,1) = 2*(c*d+a*b)   ;  
  R(0,2) = 2*(b*d+a*c)  ; R(1,2) = 2*(c*d-a*b)     ; R(2,2) = 1-2*b*b-2*c*c ;
  
  return( R );
  
}

// Function transforming a rotation matrix into a quaternion
// [[Rcpp::export]]
arma::vec Rot2QuatC (arma::mat R) {
  arma::vec q(4);
  arma::vec divisor(4);
  int m;
  
  if( R(0,0)+R(1,1)+R(2,2)==3 ){
    q(0) = 1; q(1) = 0, q(2) = 0; q(3) = 0;
  }else{  
    divisor(0) = sqrt( 1+R(0,0)+R(1,1)+R(2,2) )/2;
    divisor(1) = sqrt( 1+R(0,0)-R(1,1)-R(2,2) )/2;
    divisor(2) = sqrt( 1-R(0,0)+R(1,1)-R(2,2) )/2;
    divisor(3) = sqrt( 1-R(0,0)-R(1,1)+R(2,2) )/2;
    m = divisor.index_max();

    if( m == 0 ) {
      q(0) = divisor(m);
      q(1) = ( R(2,1) - R(1,2) ) / ( 4*divisor(m) );
      q(2) = ( R(0,2) - R(2,0) ) / ( 4*divisor(m) );
      q(3) = ( R(1,0) - R(0,1) ) / ( 4*divisor(m) );
    } else if ( m == 1){
      q(0) = ( R(2,1) - R(1,2) ) / ( 4*divisor(m) );
      q(1) = divisor(m);
      q(2) = ( R(0,1) + R(1,0) ) / ( 4*divisor(m) );
      q(3) = ( R(0,2) + R(2,0) ) / ( 4*divisor(m) );
    } else if ( m == 2 ){
      q(0) = ( R(0,2) - R(2,0) ) / ( 4*divisor(m) );
      q(1) = ( R(0,1) + R(1,0) ) / ( 4*divisor(m) );
      q(2) = divisor(m);
      q(3) = ( R(1,2) + R(2,1) ) / ( 4*divisor(m) );
    } else {
      q(0) = ( R(1,0) - R(0,1) ) / ( 4*divisor(m) );
      q(1) = ( R(0,2) + R(2,0) ) / ( 4*divisor(m) );
      q(2) = ( R(1,2) + R(2,1) ) / ( 4*divisor(m) );
      q(3) = divisor(m);
    }
  }
  if( q(0)<0 ){
    q=-q;
  }
  return( q / sqrt(sum(q%q)) );
}

////////////////////////////////////////////////////////////////
//            Geometric Operations on SO(3) and so(3)
////////////////////////////////////////////////////////////////
// Function computing Intrinsic distance between rotations
// [[Rcpp::export]]
double IntrinsicDistSO3C( arma::vec q, arma::vec p){
  double phi;
  if( sum(q%p)>=1 ){
    phi = 0;
  } else if( sum(q%p)<= -1){
    phi = M_PI;
  }else{
    phi = acos(sum(q%p));
  }
  return( 2 * min(phi, M_PI-phi) );
}

//Function computing Exponential so(3) -> SO(3)
// [[Rcpp::export]]
arma::mat ExpSO3C( arma::mat A ) {
  double absA;
  arma::mat R;
  
  if( sqrt(sum(diagvec(A.t() * A)) / 2) != 0 ){
      absA = sqrt(sum(diagvec(A.t() * A)) / 2);
      R    = eye(3,3) + sin(absA) / absA * A + (1-cos(absA)) / pow(absA,2) * A * A;
  }else{
      R = eye(3,3);
  }

  return ( R );
}

//Function computing Logarithm SO(3) -> so(3)
// [[Rcpp::export]]
arma::mat LogSO3C (arma::mat R){
  double theta = (sum(diagvec(R))-1)/2;
  arma::mat X;
  
  if( theta >= 1 & theta < 1 + 1e-4){
    X = eye(3,3);
  }else{
    double beta = acos(theta);
    X = beta / 2 / sin(beta) * (R-R.t());
  }
  
  return ( X );
}
// Function transforming vector into square matrix
// [[Rcpp::export]]
arma::mat Vec2screwMC (double xx, double yy, double zz){
  arma::mat x(3,3);
  
  x(0,0) =  0;
  x(1,0) =  zz;
  x(2,0) = -yy;
  x(0,1) = -zz;
  x(1,1) =  0;
  x(2,1) =  xx;
  x(0,2) =  yy;
  x(1,2) = -xx;
  x(2,2) =  0;
  
  return( x  );
  
}
// Function transforming screw matrix into vector
// [[Rcpp::export]]
arma::mat ScrewM2VecC (arma::mat A){
  arma::vec x(3);
  
  x(0) =  A(2,1);
  x(1) =  A(0,2);
  x(2) =  A(1,0);
  
  return( x  );
}

// Function computing the length loss between two functions sampled on the same grid
// [[Rcpp::export]]
double distFT1C ( arma::mat f, arma::mat g ) {
  
  int N = f.n_cols ;
  arma::vec p1(4) ;
  arma::vec p2(4) ;
  
  double LOSS = 0;
    
  for(int n=0; n < N-1; n++){
    p1 = QuatMultC2( f.col(n),  QuatInvC2(g.col(n) ) );
    p2 = QuatMultC2( f.col(n+1),  QuatInvC2(g.col(n+1)) ) ;

    p1 = p1 / sqrt( sum(p1%p1) )  ;
    p2 = p2 / sqrt( sum(p2%p2) )  ;
    
    LOSS = LOSS + IntrinsicDistSO3C(p1,p2) ;
  }

  return(LOSS) ;
}

// Function computing the length loss between two functions sampled on the same grid
// [[Rcpp::export]]
double distFT2C ( arma::mat f, arma::mat g ) {
  
  int N = f.n_cols ;
  arma::vec p1(4) ;
  arma::vec p2(4) ;
  
  double LOSS = 0;
  
  for(int n=0; n < N-1; n++){
    p1 = QuatMultC2( QuatInvC2(f.col(n)),  g.col(n) );
    p2 = QuatMultC2( QuatInvC2(f.col(n+1)), g.col(n+1) ) ;
    
    p1 = p1 / sqrt( sum(p1%p1) )  ;
    p2 = p2 / sqrt( sum(p2%p2) )  ;
    
    LOSS = LOSS + IntrinsicDistSO3C(p1,p2) ;
  }
  
  return(LOSS) ;
}

// Function computing the length loss between two functions sampled on the same grid
// [[Rcpp::export]]
double distFTC ( arma::mat f, arma::mat g ) {
  return( (distFT2C(f,g) + distFT1C(f,g)) / 2 );
}

////////////////////////////////////////////////////////////////
//            Ziezold Operations
////////////////////////////////////////////////////////////////
//Function computing Optimal Position of vectors
// [[Rcpp::export]]
arma::mat OptPosQuatC (arma::mat A, arma::vec qRef){
    int n = A.row(1).n_elem;
  
    if( n > 1){
      for( int i = 0; i <= n-1; i++){
        if( sum((qRef-A.col(i))%(qRef-A.col(i))) > sum((qRef+A.col(i))%(qRef+A.col(i)))) A.col(i) = -A.col(i);
      }
    }
    return( A );
}

//Function computing ProjectiveMean of a sample
// [[Rcpp::export]]
arma::vec ProjectiveMeanC (arma::mat L, int MaxIt, double err){
  arma::vec x ;
  x.zeros() ;
  bool cont   = TRUE;
  int n       = 0;
  
  if( L.n_cols > 1 ){
      arma::vec x = L.col(1);

      arma::vec xOld;
      arma::mat OptPos;
      while ( cont==TRUE ) {
        n      = n + 1;
        xOld   = x;
        OptPos = OptPosQuatC(L, x);
        x      = sum(OptPos,1);
        x      = x / norm(x);
          
            if (sqrt(sum((x-xOld)%(x-xOld))) < err || sqrt(sum((x+xOld)%(x+xOld))) < err || n > MaxIt){ cont = FALSE; }
            cont= FALSE;
      }
      if (cont==FALSE)
      { 
        return( x );
      } else {
        return( 0*ones(5) );
      }
  }else{
    return(L) ;
  }
}

// RotEstimC
// M              Number of simulations
// perm           Vector of permutations
// permC          Complement permutations of perm
// data.aligned   Aligned data for permutation test
// [[Rcpp::export]]
arma::mat RotEstimC( arma::mat A, arma::mat B ){
  arma::mat AB = A * B.t()  ;
  arma::mat u;
  arma::vec s;
  arma::mat v;

  arma::mat R = eye(4,4);
  arma::mat S = eye(4,4);

  arma::svd(u, s, v, AB);

  if( det(u)*det(v) > 0 ){
    R = u * S * v.t()  ;
  }else{
    S(3,3) = -1;
    R = u * S * v.t()  ;
  }

  return(R);
}

// GeodesicInterpolationC
// f              Matrix 4xN containing the values of the function
// times          Vector of length N containing ordered time points on which f is sampled
// new_times      Vector containing the time points on which the interpolate of f has to be evaluated
// [[Rcpp::export]]
NumericMatrix GeodesicInterpolationC( NumericMatrix f, NumericVector times, NumericVector new_times ){
 int N         = f.ncol() ;
 int M         = new_times.length() ;
 NumericMatrix inter_f(4,M) ;
 NumericVector x(4) ;
 double Normx ;


 int count = 0 ;

  for( int m = 0; m<M; m++ ){
   //find the intervall into which new_times(m) belongs

   bool found = false ;

   while( count < N && found==false ){
     if( new_times(m) >= times(count) && new_times(m) <= times(count+1)  ){
       found = true ;
     }else{
       count = count + 1 ;
     }
   }

   double dt = (new_times(m) - times(count)) / (times(count+1) - times(count));

   x = (1-dt)*f(_,count) + f(_,count+1) ;

   Normx = 0 ;
   for(int k = 0; k<4; k++){
     Normx = Normx + x(k) * x(k) ;
   }
   inter_f(_,m) = x / sqrt(Normx);
 }

  return(inter_f);
}

////////////////////////////////////////////////////////////////
//            Time Warping Dynamical Program
////////////////////////////////////////////////////////////////
//Function computing the loss internally
// [[Rcpp::export]]s
double wFunC ( NumericMatrix f, int k, NumericMatrix g, int n, int m ) {
  NumericVector w1 = g.column(n) ; NumericVector w2 = g.column(m) ;
  NumericVector v1 = f.column(k) ; NumericVector v2 = f.column(k+1) ;
  
  NumericVector p1 = QuatMultC(w1, QuatInvC(v1)) ;
  NumericVector p2 = QuatMultC(w2, QuatInvC(v2)) ;
  p1 = p1 / sqrt( sum(p1*p1) )  ;
  p2 = p2 / sqrt( sum(p2*p2) )  ;
  
  NumericVector q1 = QuatMultC( QuatInvC(w1), v1) ;
  NumericVector q2 = QuatMultC( QuatInvC(w2), v2) ;
  q1 = q1 / sqrt( sum(q1*q1) )  ;
  q2 = q2 / sqrt( sum(q2*q2) )  ;
  
  return( (IntrinsicDistSO3C( q1, q2 ) + IntrinsicDistSO3C( p1, p2 )) / 2 ) ;
}

//Function computing the dynamical program matrix
double timeWeight ( NumericVector times, int n, int m, double R ) {
  double dt = R ;
  double weight = 0 ;
  
  if ( abs( times( m ) - times( n ) ) > R ){
    weight = pow( abs( times( m ) - times( n ) ) - R, 1.1 ) ;
  }

  return( weight ) ;
}
//Function computing the dynamical program matrix
// [[Rcpp::export]]
List optimV ( NumericMatrix f, NumericMatrix g, NumericMatrix V, NumericVector times_g, double b, double R ) {
  // Maximal indices for row and column of Dijkstra-Algorithm matrix
  int N = f.ncol() - 1 ;
  int M = g.ncol() - 1 ;

  NumericMatrix Vopt = V ;
  NumericMatrix wMin(N + 1,M + 1) ;
  
  // Computing the Dijkstra-Algorithm matrix for the inner points
  for( int n = 0; n <= N - 1; n++ ){
    for( int m = n; m <= M - N + n;  m++ ){
      if (n==0) {
        Vopt( n + 1, m + 1 ) = wFunC( f, n, g, 0, m + 1 ) + b * timeWeight( times_g, 0, m, R ) ;
      }else {
        NumericVector v( m - n + 1 ) ;
        for( int l = n; l <=  m;  l++ ){
            v( l - n ) = Vopt( n, l ) + wFunC( f, n, g, l, m + 1 )
                                      + b * timeWeight( times_g, l, m, R ) ;
        }
        Vopt( n + 1, m + 1 )   = min( v ) ;
        wMin( n + 1, m + 1 )   = which_min( v ) + n ;
      }
    }
  }
  
//  for(int m = N - 1; m <= M - 1; m++ ){
//    Vopt( N, m + 1 ) = Vopt( N - 1, m ) + wFunC( f, N - 1, g, m , M  )
//                                        + b * timeWeight( times_g, m, M - 1, R ) ;
//  }
  
  List ret ;
  ret["V"] = Vopt ;
  ret["wMin"] = wMin ;
  return(ret) ;
}


//Function computing the time warped versions of two functions defined on a grid with odd(!) number of points Moreover N==M required
// [[Rcpp::export]]
NumericMatrix timeWarp ( NumericMatrix f, NumericMatrix g, double b ) {
  int N         = f.ncol() ;
  int M         = g.ncol() ;
  int reduced_N = (N+1)/2 ;

  NumericMatrix reduced_f(4,reduced_N) ;
  NumericMatrix new_g(4,N) ;
  
  // Create gGrid Vector
  NumericVector gGrid(M) ;
  gGrid(0) = 0;
  for( int m = 1; m < M; m++){
    gGrid(m) = (double) m / ((double) M-1);
  }
  
  int count = 0 ;
  for( int n=0; n < N; n=n+2 ){
    reduced_f(_,count) = f(_,n);
    count++ ;
  }

  // Initialize matrix V.
  NumericMatrix V(reduced_N,M) ;
  V.fill(R_PosInf) ;
  NumericVector v(M) ;
  v.fill(0) ;
  V(0,_) = v ;
  double R = 1 / ((double) reduced_N - 1) ; 
  List W = optimV( reduced_f, g, V, gGrid, b = b, R=R ) ;
  NumericMatrix wMin  = W["wMin"] ;
  NumericMatrix V2  = W["V"] ;
  
  // calculate gamma 
  NumericVector g_index(reduced_N) ;
  g_index.fill(0) ;
  g_index(0) = 0 ;
  g_index(reduced_N-1) = M-1 ;
  g_index(reduced_N-2) = wMin(reduced_N-1,M-1) ;
  
  for( int n = reduced_N-3; n>0; n = n-1 ){
    g_index( n ) = wMin( n+1, g_index(n+1) ) ;
  }
  
  for(int n = 1; n<N; n=n+2 ){
    new_g(_,n-1) = g(_,g_index( (n-1)/2 )) ;
    new_g(_,n+1) = g(_,g_index( (n+1)/2 )) ;
    
    NumericVector inter_g = new_g(_,n-1) + new_g(_,n+1) ;
    double NormInter_g = 0 ;
    for(int k = 0; k<4; k++){
      NormInter_g = NormInter_g + inter_g(k) * inter_g(k) ;
    }
    new_g(_,n) = inter_g / sqrt(NormInter_g);
  }
  
  return(new_g) ;
}

//Function computing the time warped versions of two functions defined on the same grid
// [[Rcpp::export]]
NumericMatrix timeWarp2 ( NumericMatrix f, NumericMatrix g, int factorN2M, double b ) {
  int N         = f.ncol() ;
  int M         = N + factorN2M * (N-1) ;
  
  NumericVector times     =  armaMat2NumMat(linspace( 0, 1, N ));
  NumericVector times_new =  armaMat2NumMat(linspace( 0, 1, M ));
  
  NumericMatrix stretched_g = GeodesicInterpolationC( g, times, times_new ) ;
  NumericMatrix new_g(4,N) ;
  
  // Initialize matrix V.
  NumericMatrix V(N,M) ;
  V.fill(R_PosInf) ;
  NumericVector v(M) ;
  v.fill(0) ;
  V(0,_) = v ;
  List W = optimV( f, stretched_g, V, times_new, b, times(1) ) ;
  NumericMatrix wMin  = W["wMin"] ;
  NumericMatrix V2  = W["V"] ;
  
  // calculate gamma
  NumericVector g_index(N) ;
  g_index.fill(0) ;
  g_index(0) = 0 ;
  g_index(N-1) = M-1 ;
  g_index(N-2) = wMin(N-1,M-1) ;

  for( int n = N-3; n>0; n = n-1 ){
    g_index( n ) = wMin( n+1, g_index(n+1) ) ;
  }

  for( int n = 0; n < N; n++ ){
    new_g(_,n) = stretched_g(_,g_index( n )) ;
  }
  
  return(new_g) ;
}

////////////////////////////////////////////////////////////////
//           Functions for Exponential Model
////////////////////////////////////////////////////////////////
// Function computing the residual coordinates for faster simulation
// [[Rcpp::export]]
List rSampleExp (NumericMatrix x, NumericMatrix y, NumericMatrix z){
  int T     = x.ncol();
  int nSamp = x.nrow();
  arma::mat mean(4,T);
  
  cube rSample(nSamp, 4, T);
  cube Residuum(nSamp, 3, T);

  NumericVector quat(4);
  
  for( int n = 0; n <= nSamp-1; n++ ){
    for( int t = 0; t <= T-1; t++ ){
      quat = Euler2QuatC(Rot2EulerC( ExpSO3C( Vec2screwMC(x(n,t),y(n,t),z(n,t)) ) ));
        for( int i = 0; i < 4; i++){
          rSample.slice(t)(n,i)  = quat(i);
        }
    }
  }
  
  for( int t = 0; t < T; t++ ){
    mean.col(t) = ProjectiveMeanC( rSample.slice(t).t(), 100, 10^(-10) );
  }
  
  NumericVector val(3);
  arma::vec q(4);
  
  for(int n = 0; n < nSamp; n++){
    for( int t = 0; t < T; t++){
      q(0) = rSample.subcube(n,0,t,n,3,t)(0,0,0);
      q(1) = rSample.subcube(n,0,t,n,3,t)(0,1,0);
      q(2) = rSample.subcube(n,0,t,n,3,t)(0,2,0);
      q(3) = rSample.subcube(n,0,t,n,3,t)(0,3,0);
      
      val = ScrewM2VecC(LogSO3C(Quat2RotC( QuatMultC2(QuatInvC2(mean.col(t)), q ) )));
      for( int i = 0; i<3; i++ ){
        Residuum.slice(t)(n,i) = val(i,0);
      }
    }
  }
  
  List ret ;
  ret["Residuum"] = Residuum ;
  ret["Mean"] = mean ;
  return(ret) ;
}

////////////////////////////////////////////////////////////////
//           Function for Simulation of covering rates
////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::cube array2cube( SEXP myArray ) {
  
  Rcpp::NumericVector vecArray(myArray);
  Rcpp::IntegerVector arrayDims = vecArray.attr("dim");
  
  arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  
  return(cubeArray);  
  
}
// Loop in permutation test
// M              Number of simulations
// perm           Vector of permutations
// permC          Complement permutations of perm
// data.aligned   Aligned data for permutation test
// [[Rcpp::export]]
NumericVector PermTestLoop ( int M, arma::cube dataAligned, IntegerMatrix perm, IntegerMatrix permC ) {
      int T  = dataAligned.n_cols ;
      int N1 =  perm.nrow() ;
      int N2 =  permC.nrow() ;
      
      arma::cube Samp1(4, N1, T ) ;
      arma::cube Samp2(4, N2, T ) ;
      
      arma::mat mean1(4,T) ;
      arma::mat mean2(4,T) ;
      
      NumericVector distVec(M) ;
      
      int N = dataAligned.n_slices ;
      int count = 0 ;
      
    for( int m = 0; m < M; m++ ){
        count = 0 ;
        for(int k = 0; k < N; k++){
          if( k == perm(count, m)-1 ){
              Samp1.tube(span::all, span(count,count)) = dataAligned.slice(k) ;
              count = count + 1 ;
          }
        }
        
        count = 0 ;
        for(int k = 0; k < N; k++){
          if( k == permC(count, m)-1 ){
                Samp2.tube(span::all, span(count,count)) = dataAligned.slice(k) ;
                count = count + 1 ;
          }
        }
        
        for( int t = 0; t < T; t++ ){
           mean1.col(t) = ProjectiveMeanC( Samp1.slice(t), 100, 10^(-10) ) ;
           mean2.col(t) = ProjectiveMeanC( Samp2.slice(t), 100, 10^(-10) ) ;
        }
     
        // NumericVector times     =  armaMat2NumMat(linspace( 0, 1, T ));
        // NumericVector times_new =  armaMat2NumMat(linspace( 0, 1, 3*T-2 ));
        // distVec(m) = distFTC(NumMat2armaMat(GeodesicInterpolationC( armaMat2NumMat(mean1), times, times_new)), NumMat2armaMat(GeodesicInterpolationC( armaMat2NumMat(mean2), times, times_new))) ;
     
       distVec(m) = distFTC(mean1, mean2) ;
    }
    
    return(distVec) ;
}

// Loop in permutation test
// M              Number of simulations
// perm           Vector of permutations
// permC          Complement permutations of perm
// data.aligned   Aligned data for permutation test
// [[Rcpp::export]]
arma::mat PermTest2Loop ( int M, arma::cube dataAligned, arma::mat perm, arma::mat permC ) {
  int T  = dataAligned.n_cols ;
  int N1 =  perm.n_rows ;
  int N2 =  permC.n_rows ;

  int N11 = 0 ;
  int N21 = 0 ;

  arma::mat mean11(4,T) ;
  arma::mat mean12(4,T) ;
  arma::mat mean21(4,T) ;
  arma::mat mean22(4,T) ;

  arma::mat R1 = eye(4,4) ;
  arma::mat R2 = eye(4,4) ;

  arma::mat mean1(4,T) ;
  arma::mat mean2(4,T) ;

  NumericVector distVec(M) ;

  int N = dataAligned.n_slices ;
  int count = 0 ;

  for( int m = 0; m < M; m++ ){
        N11 = 0 ;
        N21 = 0 ;
    
        for(int l = 0; l < N1; l++){
          if( perm(l,m) <= N1 ){
            N11 = N11 + 1 ;
          }
        }
        for(int l = 0; l < N2; l++){
          if( permC(l,m) <= N1 ){
            N21 = N21 + 1 ;
          }
        }
        
        arma::cube Samp11(4, N11, T )       ;
        arma::cube Samp12(4, N1 - N11, T )  ;
        arma::cube Samp21(4, N21, T )       ;
        arma::cube Samp22(4, N2 - N21, T )  ;
    
        Samp11.zeros() ;
        Samp12.zeros() ;
        Samp21.zeros() ;
        Samp22.zeros() ;
    
        mean11.zeros() ;
        mean12.zeros() ;
        mean21.zeros() ;
        mean22.zeros() ;
    
        count = 0 ;
        for( int k = 0; k < N1; k++ ){
          if( perm(k, m) <= N1 ){
            Samp11.tube(span::all, span(count,count)) = dataAligned.slice( perm(k, m) - 1 ) ;
            count = count + 1 ;
          }
        }
        count = 0 ;
        for(int k = 0; k < N1; k++){
          if( perm(k, m) > N1 ){
            Samp12.tube(span::all, span(count,count)) = dataAligned.slice( perm(k, m) - 1 ) ;
            count = count + 1 ;
          }
        }
    
        count = 0 ;
        for(int k = N1; k < N; k++){
          if( permC(k-N1, m) <= N1 ){
            Samp21.tube(span::all, span(count,count)) = dataAligned.slice( permC(k-N1, m) - 1 ) ;
            count = count + 1 ;
          }
        }
        count = 0 ;
        for(int k = N1; k < N; k++){
          if( permC(k-N1, m) > N1 ){
            Samp22.tube(span::all, span(count,count)) = dataAligned.slice( permC(k-N1, m) - 1 ) ;
            count = count + 1 ;
          }
        }
    
        for( int t = 0; t < T; t++ ){
          if( N11 != 0 ){
              mean11.col(t) = ProjectiveMeanC( Samp11.slice(t), 100, 10^(-10) ) ;
          }
          if( N1 - N11 != 0 ){
            mean12.col(t) = ProjectiveMeanC( Samp12.slice(t), 100, 10^(-10) ) ;
          }
          if( N21 != 0 ){
            mean21.col(t) = ProjectiveMeanC( Samp21.slice(t), 100, 10^(-10) ) ;
          }
          if( N2 - N21 != 0 ){
            mean22.col(t) = ProjectiveMeanC( Samp22.slice(t), 100, 10^(-10) ) ;
          }
        }
    
        // catch the case that one sample consists entirely of one group
        if( N11 != 0 && N1 - N11 != 0 ){
          R1 = RotEstimC( mean11, mean12 ) ;
        }else{
          R1 = eye(4,4) ;
        }
        if( N21 != 0 &&  N2 - N21 != 0 ){
          R2 = RotEstimC( mean21, mean22 );
        }else{
          R2 = eye(4,4) ;
        }
    
        for(int t = 0; t < T; t++ ){
          mean12.col(t) = R1 * mean12.col(t) ;
          mean22.col(t) = R2 * mean22.col(t) ;
    
          arma::vec x ;
          if( norm(mean12.col(t) - mean11.col(t)) <= norm(mean12.col(t) + mean11.col(t)) ){
            x = (mean12.col(t) + mean11.col(t)) / 2 ;
          }else{
            x = (mean12.col(t) - mean11.col(t)) / 2 ;
          }
          mean1.col(t) = x / norm(x) ;
    
          if( norm(mean22.col(t) - mean21.col(t)) <= norm(mean22.col(t) + mean21.col(t)) ){
            x = (mean22.col(t) + mean21.col(t)) / 2 ;
          }else{
            x = (mean22.col(t) - mean21.col(t)) / 2 ;
          }
          mean2.col(t) = x / norm(x) ;
        }
    
        arma::mat R = RotEstimC( mean1, mean2 );
        for(int t = 0; t < T; t++ ){
          mean2.col(t) = R * mean2.col(t) ;
        }
    
        distVec(m) = distFTC(mean1, mean2) ;
  }
  
  return( distVec ) ;
}


// Loop in permutation test 3
// M              Number of simulations
// perm           Vector of permutations
// permC          Complement permutations of perm
// data.aligned   Aligned data for permutation test
// [[Rcpp::export]]
arma::vec PermTest3Loop ( int M, arma::cube dataAligned, arma::mat perm, arma::mat permC, int N2Mfac, double b ) {
  int T  = dataAligned.n_cols ;
  int N1 =  perm.n_rows ;
  int N2 =  permC.n_rows ;
  
  int N11 = 0 ;
  int N21 = 0 ;
  
  arma::mat mean11(4,T) ;
  arma::mat mean12(4,T) ;
  arma::mat mean21(4,T) ;
  arma::mat mean22(4,T) ;
  
  arma::mat R1 = eye(4,4) ;
  arma::mat R2 = eye(4,4) ;
  
  arma::mat mean1(4,T) ;
  arma::mat mean2(4,T) ;
  
  NumericVector distVec(M) ;
  
  int N = dataAligned.n_slices ;
  int count = 0 ;
  
  for( int m = 0; m < M; m++ ){
    N11 = 0 ;
    N21 = 0 ;
    
    for(int l = 0; l < N1; l++){
      if( perm(l,m) <= N1 ){
        N11 = N11 + 1 ;
      }
    }
    for(int l = 0; l < N2; l++){
      if( permC(l,m) <= N1 ){
        N21 = N21 + 1 ;
      }
    }
    
    arma::cube Samp11(4, N11, T )       ;
    arma::cube Samp12(4, N1 - N11, T )  ;
    arma::cube Samp21(4, N21, T )       ;
    arma::cube Samp22(4, N2 - N21, T )  ;
    
    Samp11.zeros() ;
    Samp12.zeros() ;
    Samp21.zeros() ;
    Samp22.zeros() ;
    
    mean11.zeros() ;
    mean12.zeros() ;
    mean21.zeros() ;
    mean22.zeros() ;
    
    count = 0 ;
    for( int k = 0; k < N1; k++ ){
      if( perm(k, m) <= N1 ){
        Samp11.tube(span::all, span(count,count)) = dataAligned.slice( perm(k, m) - 1 ) ;
        count = count + 1 ;
      }
    }
    count = 0 ;
    for(int k = 0; k < N1; k++){
      if( perm(k, m) > N1 ){
        Samp12.tube(span::all, span(count,count)) = dataAligned.slice( perm(k, m) - 1 ) ;
        count = count + 1 ;
      }
    }
    
    count = 0 ;
    for(int k = N1; k < N; k++){
      if( permC(k-N1, m) <= N1 ){
        Samp21.tube(span::all, span(count,count)) = dataAligned.slice( permC(k-N1, m) - 1 ) ;
        count = count + 1 ;
      }
    }
    count = 0 ;
    for(int k = N1; k < N; k++){
      if( permC(k-N1, m) > N1 ){
        Samp22.tube(span::all, span(count,count)) = dataAligned.slice( permC(k-N1, m) - 1 ) ;
        count = count + 1 ;
      }
    }
    
    for( int t = 0; t < T; t++ ){
      if( N11 != 0 ){
        mean11.col(t) = ProjectiveMeanC( Samp11.slice(t), 100, 10^(-10) ) ;
      }
      if( N1 - N11 != 0 ){
        mean12.col(t) = ProjectiveMeanC( Samp12.slice(t), 100, 10^(-10) ) ;
      }
      if( N21 != 0 ){
        mean21.col(t) = ProjectiveMeanC( Samp21.slice(t), 100, 10^(-10) ) ;
      }
      if( N2 - N21 != 0 ){
        mean22.col(t) = ProjectiveMeanC( Samp22.slice(t), 100, 10^(-10) ) ;
      }
    }
    
    // catch the case that one sample consists entirely of one group
    if( N11 != 0 && N1 - N11 != 0 ){
      R1 = RotEstimC( mean11, mean12 ) ;
    }else{
      R1 = eye(4,4) ;
    }
    if( N21 != 0 &&  N2 - N21 != 0 ){
      R2 = RotEstimC( mean21, mean22 );
    }else{
      R2 = eye(4,4) ;
    }
    
    for(int t = 0; t < T; t++ ){
      mean12.col(t) = R1 * mean12.col(t) ;
      mean22.col(t) = R2 * mean22.col(t) ;
    }
    
    if( N11 != 0 && N1 - N11 != 0 ){
      mean12 = NumMat2armaMat( timeWarp2 ( armaMat2NumMat(mean11), armaMat2NumMat(mean12), N2Mfac, b ) );
    }
    if( N21 != 0 &&  N2 - N21 != 0 ){
      mean22 = NumMat2armaMat( timeWarp2 ( armaMat2NumMat(mean21), armaMat2NumMat(mean22), N2Mfac, b ) );    
    }
    
    for(int t = 0; t < T; t++ ){            
      arma::vec x ;
      if( norm(mean12.col(t) - mean11.col(t)) <= norm(mean12.col(t) + mean11.col(t)) ){
        x = (mean12.col(t) + mean11.col(t)) / 2 ;
      }else{
        x = (mean12.col(t) - mean11.col(t)) / 2 ;
      }
      mean1.col(t) = x / norm(x) ;
      
      if( norm(mean22.col(t) - mean21.col(t)) <= norm(mean22.col(t) + mean21.col(t)) ){
        x = (mean22.col(t) + mean21.col(t)) / 2 ;
      }else{
        x = (mean22.col(t) - mean21.col(t)) / 2 ;
      }
      mean2.col(t) = x / norm(x) ;
    }
    
    arma::mat R = RotEstimC( mean1, mean2 );
    for(int t = 0; t < T; t++ ){
      mean2.col(t) = R * mean2.col(t) ;
    }
    
    mean2 = NumMat2armaMat( timeWarp2( armaMat2NumMat(mean1), armaMat2NumMat(mean2), N2Mfac, b ) );
    R = RotEstimC( mean1, mean2 );
    for(int t = 0; t < T; t++ ){
      mean2.col(t) = R * mean2.col(t) ;
    }
    
    distVec(m) = distFTC(mean1, mean2) ;
  }
  
  return( distVec ) ;
}

//Function for testing C
// [[Rcpp::export]]
arma::vec TestC ( int a, int b, int N) {
  return( linspace( a, b, N ) );
}

//Function for testing C
// [[Rcpp::export]]
int Test ( arma::vec v ) {
  return( v.index_max() );
}
