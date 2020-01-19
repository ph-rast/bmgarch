// BEKK-Parameterization
functions { 
  //#include /functions/cov2cor.stan
matrix cov2cor(matrix C){
  int dm = rows(C);
  matrix[dm,dm] s;
  matrix[dm,dm] R;
  s = diag_matrix( 1 ./ sqrt(diagonal(C)) );
  R = s*C*s ; //quad_form_diag(C, s);
  return R;
}
  
}

data {
  //#include /data/data.stan
  int<lower=1> T; // lengh of time series
  int<lower=1> nt;    // number of time series
  int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
  int<lower=1> P; // AR component in MGARCH(P,Q), matrix B
  vector[nt] rts[T];  // multivariate time-series
  vector[nt] xH[T];  // time-varying predictor for conditional H
  int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
  int<lower=0, upper=1> meanstructure; // Select model for location

}

transformed data {
  //#include /transformed_data/xh_marker.stan
  cov_matrix[nt] H1_init = diag_matrix( rep_vector(.5, nt));
}

parameters { 
  // ARMA parameters
  vector[nt] phi0; 
  corr_matrix[nt] C_R; // Const is symmetric, A, B, are not
  vector<lower= 0>[nt] C_sd;
  
  //  matrix[nt, nt] A;
  real<lower = 0, upper = 1> a_11; // one element of A needs to be non-negative: Cf Engle Kroner (1995) p. 128
  real a_12;
  real a_21;
  real<upper = 1 > a_22;
  //  matrix[nt, nt] B;
  real<lower = 0, upper = (1 - a_11) > b_11; // one element of B needs to be non-negative: Cf Engle Kroner (1995) p. 128
  real b_12;
  real b_21;
  real<upper = (1 - a_22) > b_22;
  
}

transformed parameters {
  matrix[nt, nt] A[Q];
  matrix[nt, nt] B[P];
  cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];

  matrix[nt, nt] A_part = diag_matrix( rep_vector(0.0, nt));
  matrix[nt, nt] B_part = diag_matrix( rep_vector(0.0, nt));

  cov_matrix[nt] Cnst = quad_form( diag_matrix( C_sd ), C_R );
  // redefine A and B
  A[1,1,1] = a_11;
  A[1,1,2] = a_12;
  A[1,2,1] = a_21;
  A[1,2,2] = a_22;

  B[1,1,1] = b_11;
  B[1,1,2] = b_12;
  B[1,2,1] = b_21;
  B[1,2,2] = b_22;
  
  
  
  // Initialize model parameters
  mu[1,] = phi0;
  H[1,] = H1_init;
  L_H[1,] = cholesky_decompose(H[1,]); // cf. p 69 in stan manual for how to index

  for (t in 2:T){    
    // Meanstructure model:
    //#include /model_components/mu.stan
    mu[t, ] = phi0;

    // reset A_part and B_part to zero for each iteration t
     A_part = diag_matrix( rep_vector(0.0, nt));
     B_part = diag_matrix( rep_vector(0.0, nt));
        
    for (q in 1:min( t-1, Q) ) {
      rr[t-q,] = ( rts[t-q,] - mu[t-q,] )*( rts[t-q,] - mu[t-q,] )';
       A_part = A_part + A[q]' * rr[t-q,] * A[q];
       }
      for (p in 1:min( t-1, P) ) {
       B_part = B_part + B[p]' * H[t-p,] * B[p];
        }
    //    if( xH_marker == 0 ) {
      //      C_sd = exp( beta0 ); 
      //Cnst =  quad_form_diag(C_R, C_sd );
           H[t,] = Cnst + A_part +  B_part;
	   //H[t,] = Cnst + A[1]' * rr[t-q,] * A[q] +  B[p]' * H[t-p,] * B[p];
      // } else if( xH_marker >= 1) {
      // C_sd = exp( append_col( 1.0, xH[t]' ) * beta ); 
      //Cnst =  quad_form_diag(C_R, C_sd );
      //  H[t,] = Cnst  + A_part +  B_part;
      
    L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  // https://mc-stan.org/documentation/case-studies/mle-params.html
// priors
  a_11 ~ normal(0.5, .2);
  a_12 ~ normal(0, 1);
  a_11 ~ normal(0, 1);
  a_22 ~ normal(0.5, .1);
  b_11 ~ normal(0.5, .2);
  b_12 ~ normal(0, 1);
  b_11 ~ normal(0, 1);
  b_22 ~ normal(0.5, .1);
  // Prior on nu for student_t
  //if ( distribution == 1 )
  //  nu ~ normal( nt, 50 );
  // Prior for initial state
  //H1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  // to_vector(theta) ~ normal(0, 1);
  //  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  //Cnst ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  /* to_vector(beta0) ~ normal(-2, 5); */
  /* to_vector(beta1) ~ normal(0, 1); */
  C_sd ~ normal( 1, 1);
   C_R ~ lkj_corr( 2 ); 
  

  
  // for( k in 1:nt){ 
  // for( q in 1:Q ) {
       //  A_diag[q,k] ~ uniform(0, 1);
       //target += log( A_diag[q,k] ) + log1m( A_diag[q,k] );       
       //  target += uniform_lpdf(A_diag[q,k] | 0, 1 ) +  log( A_diag[q,k] ) + log1m( A_diag[q,k] );
     /*   target += uniform_lpdf(Av[q,k] | -Ca[q,k], Ca[q,k]) + log( 2*Ca[q,k] ) + log_inv_logit( A_log[q,k] ) + log1m_inv_logit( A_log[q,k] ); */
  // }
  // for ( p in 1:P ) {
       //B_diag[p,k] ~ uniform(0, 1);
       //target += log( B_diag[p,k] ) + log1m( B_diag[p,k] );
       //target += uniform_lpdf(B_diag[p,k] | 0, 1 ) +  log( B_diag[p,k] ) + log1m( B_diag[p,k] );
       //target += uniform_lpdf(Bv[p,k] | -Cb[p,k], Cb[p,k]) + log( 2*Cb[p,k] ) + log_inv_logit( B_log[p,k] ) + log1m_inv_logit( B_log[p,k] );
  // }
  //}
  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      //      rts[t,] ~ multi_normal_cholesky(mu[t,], L_H[t,]);
      target += multi_normal_cholesky_lpdf( rts[t, ] | mu[t,], L_H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      // rts[t,] ~ multi_student_t(nu, mu[t,], L_H[t,]*L_H[t,]');
      target += multi_student_t_lpdf( rts[t, ] | 10, mu[t,], L_H[t,]*L_H[t,]');
    }
  }
}
//
generated quantities {
/*   matrix[nt,T] rts_out; */
/*   real log_lik[T]; */
/*   corr_matrix[nt] corC; */
/*   corr_matrix[nt] corH[T]; */
/* /\*   //  row_vector[nt] C_var; *\/ */

/* /\* //Const = multiply_lower_tri_self_transpose(Cnst); *\/ */
/*     corC = cov2cor(Cnst);  */
/* /\*   // Square SD to get var *\/ */
/* /\*   //  C_var = exp( beta0 ) .* exp( beta0 ); *\/ */

/* /\*   // retrodict *\/ */
/* #include /generated/retrodict_H.stan  */

}
