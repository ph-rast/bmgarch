data {
#include /data/data.stan  
}

transformed data {
#include /transformed_data/xh_marker.stan  
}


parameters {
  // ARMA parameters
#include /parameters/arma.stan
  // predictor for H
  //cov_matrix[ xH_marker >= 1 ? nt : 0 ] beta;
  cov_matrix[nt] beta;
  
  // DF constant nu for student t
  real< lower = 2 > nu;

  //  cholesky_factor_cov[nt] Cnst; // Const is symmetric, A, B, are not
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  

  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];
  matrix[nt, nt] A[Q];
  matrix[nt, nt] B[P];
}

generated quantities {
  // Define matrix for rts prediction
  vector[nt] rts_p[ahead + max(Q,P)];
  cov_matrix[nt] H_p[ahead + max(Q,P)];

  matrix[nt,nt] rr_p[ahead + max(Q,P)];
  vector[nt] mu_p[ahead + max(Q,P)];

  // Placeholders
  matrix[nt, nt] A_part_p;
  matrix[nt, nt] B_part_p;
  
  // Populate with non-NA values to avoid Error in stan
  rts_p[1:(ahead + max(Q,P)), ] = rts[ 1:(ahead + max(Q,P)), ];
  H_p[  1:(ahead + max(Q,P)), ] = H[  1:(ahead + max(Q,P)), ];
  mu_p[ 1:(ahead + max(Q,P)), ] = mu[ 1:(ahead + max(Q,P)), ];
  rr_p[ 1:(ahead + max(Q,P)), ] = rr[ 1:(ahead + max(Q,P)), ];
  
  // Obtain needed elements from mu and fill into mu_p
  rts_p[1:max(Q, P), ] = rts[ (T-(max(Q,P)-1) ):T, ];
  H_p[  1:max(Q, P), ] =  H[ (T-(max(Q,P)-1) ):T, ];
  mu_p[ 1:max(Q, P), ] = mu[ (T-(max(Q,P)-1) ):T, ];
  // rr is of length T-1
  rr_p[ 1:max(Q, P), ] = rr[ (T-1-(max(Q,P)-1) ):(T-1), ];
  
  // Forecast
  for (t in (max(Q, P) + 1 ):( max(Q, P) + ahead ) ){

    // reset both matrices to zero for each iteration t
    A_part_p = diag_matrix( rep_vector(0.0, nt));
    B_part_p = diag_matrix( rep_vector(0.0, nt));
    
    if( meanstructure == 0 ){
      mu_p[t, ] = phi0;
    } else if( meanstructure == 1 ){
      mu_p[t, ] = phi0 + phi * rts_p[t-1, ] + theta * (rts_p[t-1, ] - mu_p[t-1,]) ;
    }

    for (q in 1:min( t-1, Q) ) {
      rr_p[t-q,] = ( rts_p[t-q,] - mu_p[t-q,] )*( rts_p[t-q,] - mu_p[t-q,] )';
      A_part_p = A_part_p + A[q]' * rr_p[t-q,] * A[q];
    }
    for (p in 1:min( t-1, P) ) {
      B_part_p = B_part_p + B[p]' * H_p[t-p,] * B[p];
    }
    if( xH_marker == 0 ) {
      H_p[t,] = Cnst + A_part_p +  B_part_p;
    } else if( xH_marker >= 1) {
      H_p[t,] = Cnst + beta * xH_m[t]  + A_part_p +  B_part_p;
    } 
    
    if ( distribution == 0 ) {
      rts_p[t,] = multi_normal_rng( mu_p[t,], H_p[t,]);
    } else if ( distribution == 1 ) {
      rts_p[t,] = multi_student_t_rng( nu, mu_p[t,], H_p[t,]);
    }
  }
}
