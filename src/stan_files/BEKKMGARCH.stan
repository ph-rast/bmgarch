// BEKK-Parameterization
functions { 
#include /functions/cov2cor.stan
}

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
  row_vector[nt] beta0;
  vector[nt] beta1;
  
  // in case Cnst is predicted, separate into C_sd*C_R*C_sd
  corr_matrix[nt] C_R;

  matrix[nt, nt] A_raw[Q];
  matrix[nt, nt] B_raw[P];

    // H1 init
  cov_matrix[nt] H1_init; 
  real< lower = 2 > nu; // nu for student_t

}
transformed parameters {
  // cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];

  matrix[nt, nt] A_part = diag_matrix( rep_vector(0.0, nt));
  matrix[nt, nt] B_part = diag_matrix( rep_vector(0.0, nt));

  matrix[nt+1, nt] beta = append_row( beta0, diag_matrix(beta1) );
  row_vector[nt] C_sd;
//  cholesky_factor_cov[nt] Cnst; // Const is symmetric, A, B, are not
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  
    
  // Initialize model parameters
  mu[1,] = phi0;
  H[1,] = H1_init;
  // L_H[1,] = cholesky_decompose(H[1,]); // cf. p 69 in stan manual for how to index

  for (t in 2:T){    
    // Meanstructure model:
#include /model_components/mu.stan

    // reset A_part and B_part to zero for each iteration t
    A_part = diag_matrix( rep_vector(0.0, nt));
    B_part = diag_matrix( rep_vector(0.0, nt));
        
    for (q in 1:min( t-1, Q) ) {
      rr[t-q,] = ( rts[t-q,] - mu[t-q,] )*( rts[t-q,] - mu[t-q,] )';
      A_part = A_part + A_raw[q]' * rr[t-q,] * A_raw[q];
    }
    for (p in 1:min( t-1, P) ) {
      B_part = B_part + B_raw[p]' * H[t-p,] * B_raw[p];
    }
    if( xH_marker == 0 ) {
      C_sd = exp( beta0 ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst + A_part +  B_part;
    } else if( xH_marker >= 1) {
      C_sd = exp( append_col( 1.0, xH[t]' ) * beta ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst  + A_part +  B_part;
    }
    // L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  // https://mc-stan.org/documentation/case-studies/mle-params.html

  // Prior on nu for student_t
  nu ~ normal( nt, 50 );

  // Prior for initial state
  H1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );

  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);

  to_vector(beta0) ~ std_normal();
  to_vector(beta1) ~ std_normal();
  C_R ~ lkj_corr( 1 );

  for(q in 1:Q) {
    to_vector(A_raw[q]) ~ std_normal();
  }
  for(p in 1:P) {
    to_vector(B_raw[p]) ~ std_normal();
  }
    
  
  //for( k in 1:nt){ 
  //   for( q in 1:Q ) {
  //A_raw[q,k,k] ~ uniform(0, 1);
	//target += log( A_[q,k,k] ) + log1m( A_diag[q,k] );       
       // target += uniform_lpdf(A_diag[q,k] | 0, 1 ) +  log( A_diag[q,k] ) + log1m( A_diag[q,k] );
     /*   target += uniform_lpdf(Av[q,k] | -Ca[q,k], Ca[q,k]) + log( 2*Ca[q,k] ) + log_inv_logit( A_log[q,k] ) + log1m_inv_logit( A_log[q,k] ); */
  //  }
     // for ( p in 1:P ) {
       //B[p,k,k] ~ uniform(0, 1);
       //target += log( B[p,k,k] ) + log1m( B[p,k,k] );
       // target += uniform_lpdf(B_diag[p,k] | 0, 1 ) +  log( B_diag[p,k] ) + log1m( B_diag[p,k] );
       //target += uniform_lpdf(Bv[p,k] | -Cb[p,k], Cb[p,k]) + log( 2*Cb[p,k] ) + log_inv_logit( B_log[p,k] ) + log1m_inv_logit( B_log[p,k] );
  // }
  // }
  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      //      rts[t,] ~ multi_normal_cholesky(mu[t,], L_H[t,]);
      // target += multi_normal_cholesky_lpdf( rts[t, ] | mu[t,], L_H[t,]);
      target += multi_normal_lpdf(rts[t, ]| mu[t,], H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      // rts[t,] ~ multi_student_t(nu, mu[t,], L_H[t,]*L_H[t,]');
      target += multi_student_t_lpdf( rts[t, ] | nu, mu[t,], H[t,]);
    }
  }
}
//
generated quantities {
  matrix[nt, nt] A[Q] = A_raw;
  matrix[nt, nt] B[P] = B_raw;
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corC;
  corr_matrix[nt] corH[T];
  row_vector[nt] C_var;

  for(q in 1:Q) {
    if(A[q,1,1] < 0) {
      A[q] = -A[q];
    }
  }
  for(p in 1:P) {
    if(B[p,1,1] < 0) {
      B[p] = -B[p];
    }
  }

//Const = multiply_lower_tri_self_transpose(Cnst);
  corC = C_R;
  C_var = exp(2*beta0 );

  // retrodict
#include /generated/retrodict_H.stan

}
