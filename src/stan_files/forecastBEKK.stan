functions { 
#include /functions/cov2cor.stan
}

data {
#include /data/gq_data.stan
}

transformed data {
  matrix[ahead + max(Q,P), nt] xC_c;
#include /transformed_data/xh_marker.stan
  // Concatenate xC with xC_p 
  for(i in 1:max(Q,P) ) {
    xC_c[i] = xC[T - (max(Q,P) - 1)]';
  }
  for(i in 1:ahead) {
    xC_c[i + max(Q,P)] = xC_p[i]';
  }
}


parameters {
  // ARMA parameters
#include /parameters/arma.stan

  // predictor for H
  //cov_matrix[ xC_marker >= 1 ? nt : 0 ] beta;
  //  cov_matrix[nt] beta;
  row_vector[nt] beta0;
  vector[nt] beta1;
  
  // in case Cnst is predicted, separate into C_sd*C_R*C_sd
  corr_matrix[nt] C_R;
  
  // DF constant nu for student t
  real< lower = 2 > nu;

  //  cholesky_factor_cov[nt] Cnst; // Const is symmetric, A, B, are not
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  

  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];
  matrix[nt, nt] A[Q];
  matrix[nt, nt] B[P];
  corr_matrix[nt] corH[T];
}

generated quantities {
  // Define matrix for rts prediction
  vector[nt] rts_p[ahead + max(Q,P)];
  vector[nt] rts_forecasted[ahead];
  
  cov_matrix[nt] H_p[ahead + max(Q,P)];
  corr_matrix[nt] R_p[ahead + max(Q,P)];
  cov_matrix[nt] H_forecasted[ahead];
  corr_matrix[nt] R_forecasted[ahead];

  matrix[nt,nt] rr_p[ahead + max(Q,P)];
  vector[nt] mu_p[ahead + max(Q,P)];
  vector[nt] mu_forecasted[ahead];

  matrix[nt+1, nt] beta = append_row( beta0, diag_matrix(beta1) );

  // log lik for LFO-CV
  // only compute log_lik if it is actually requested
  real log_lik[compute_log_lik ==1 ? ahead:0];
  
  // Placeholders
  matrix[nt, nt] A_part_p;
  matrix[nt, nt] B_part_p;
  
  // Populate with non-NA values to avoid Error in stan
  rts_p[1:(ahead + max(Q,P)), ] = rts[ 1:(ahead + max(Q,P)), ];
  H_p[  1:(ahead + max(Q,P)), ] = H[  1:(ahead + max(Q,P)), ];
  mu_p[ 1:(ahead + max(Q,P)), ] = mu[ 1:(ahead + max(Q,P)), ];
  rr_p[ 1:(ahead + max(Q,P)), ] = rr[ 1:(ahead + max(Q,P)), ];
  
  R_p[ 1:(ahead + max(Q,P)), ] = corH[ 1:(ahead + max(Q,P)), ];
  
  // Obtain needed elements from mu and fill into mu_p
  rts_p[1:max(Q, P), ] = rts[ (T-( max(Q,P)-1 ) ):T, ];
  H_p[  1:max(Q, P), ] =  H[ (T-( max(Q,P)-1 ) ):T, ];
  mu_p[ 1:max(Q, P), ] = mu[ (T-( max(Q,P)-1 ) ):T, ];
  // rr is of length T-1
  rr_p[ 1:max(Q, P), ] = rr[ (T-1-( max(Q,P)-1 ) ):(T-1), ];
  R_p[  1:max(Q, P), ] =  corH[ (T - ( max(Q,P)-1 ) ):T, ];

  
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
    if( xC_marker == 0 ) {
      H_p[t,] = quad_form_diag(C_R, exp( beta0 ) ) + A_part_p +  B_part_p;
    } else if( xC_marker >= 1) {
      H_p[t,] = quad_form_diag(C_R,exp( append_col( 1.0, xC_c[ t, ] ) * beta )) +
	A_part_p +  B_part_p;
    }
    R_p[t, ] = cov2cor(H_p[t,]);

    /* sampling distributions */
#include /generated/forecast_sampling.stan
    /* if ( distribution == 0 ) { */
    /*   rts_p[t,] = multi_normal_rng( mu_p[t,], H_p[t,]); */
    /*   if( compute_log_lik ) { */
    /* 	for( i in 1:ahead ){ */
    /* 	  log_lik[i] = multi_normal_lpdf( future_rts[i] |  mu_p[t,], H_p[t,] ); */
    /* 	} */
    /*   } */
    /* } else if ( distribution == 1 ) { */
    /*   rts_p[t,] = multi_student_t_rng( nu, mu_p[t,], H_p[t,]); */
    /*   if( compute_log_lik ) { */
    /* 	for( i in 1:ahead ){ */
    /* 	  log_lik[i] = multi_student_t_lpdf( future_rts[i] | nu, mu_p[t,], H_p[t,] ); */
    /* 	} */
    /*   } */
    /* } */
  }
  rts_forecasted = rts_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  H_forecasted = H_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  R_forecasted = R_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  mu_forecasted = mu_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
#include /generated/forecast_log_lik.stan
}
