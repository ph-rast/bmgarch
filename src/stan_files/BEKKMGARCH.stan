// BEKK-Parameterization
functions { 
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
  int<lower=2> T;
  int<lower=1> nt;    // number of time series
  vector[nt] rts[T];  // multivariate time-series
  int<lower=0> ahead; // how many ahead predictions 
  matrix[nt,nt] sigma1;
}
transformed data {
  // Reverse the rts vector
  vector[nt] rev[T];
  for( i in 1:nt ) {
    rev[,i] = rts[,nt-i+1];
  }
}
parameters { 
  //  cholesky_factor_cov[nt] Cnst; // Const is symmetric, A, B, are not
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  
  // construct A, so that one element (a11) can be constrained to be non-negative, instead of adding prior onto A[1,1]
  real<lower = 0, upper = 1> Ap11;
  row_vector[nt-2] Ap1k;
  matrix[nt-1, nt-1] Ap_sub;
  //
  real<lower = 0, upper = 1> Bp11;
  row_vector[nt-2] Bp1k;
  matrix[nt-1, nt-1] Bp_sub;
  //
  vector[nt] A_log;
  vector[nt] B_log;
  // ARMA parameters
  vector[nt] phi0; 
  matrix[nt,nt] phi;
  matrix[nt,nt] theta;
}
transformed parameters {
  cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];
  matrix[nt, nt] A;
  matrix[nt, nt] B;
  vector[nt] Ca; // Upper (and lower) boundary for A 
  vector[nt] Av;
  vector[nt] Cb; // Upper (and lower) boundary for A 
  vector[nt] Bv; 
  matrix[nt, nt -1 ] Ap;
  matrix[nt, nt -1 ] Bp;
  Ap = append_row(append_col(Ap11, Ap1k), Ap_sub);
  Bp = append_row(append_col(Bp11, Bp1k), Bp_sub);
  for ( i in 1:nt) {
    Ca[i] = sqrt( 1 - dot_self(Ap[i]) );
    Av[i] = -Ca[i] + 2*Ca[i] * inv_logit( A_log[i] );
    Cb[i] = sqrt( 1 - dot_self(Bp[i]) );
    Bv[i] = -Cb[i] + 2*Cb[i] * inv_logit( B_log[i] );
  }  
  A = append_col(Ap, Av);
  B = append_col(Bp, Bv);
  // Initialize
  H[1,] = sigma1;              //rts[,1]*transpose(rts[,1]); // Initial state
  L_H[1,] = cholesky_decompose(H[1,]); // cf. p 69 in stan manual for how to index
  // = AR + MA 
  mu[1,] = phi0 + phi * rts[1, ] + (rts[1, ] - phi0) - theta * (rts[1, ] - phi0) ;
  //
  for (t in 2:T){
    mu[t, ] = phi0 + phi * rts[t-1, ] + (rts[t-1, ] - mu[t-1,]) - theta * (rts[t-1, ] - mu[t-1,]) ;
    rr[t-1,] = ( rts[t-1,] - mu[t-1,] )*( rts[t-1,] - mu[t-1,] )';
    //  H[t,] = multiply_lower_tri_self_transpose(Cnst) + A' * rr[t-1,] * A + B' * H[t-1,] * B;
    H[t,] = Cnst + A' * rr[t-1,] * A + B' * H[t-1,] * B;    
    L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  Cnst ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  for( k in 1:nt){
    target += uniform_lpdf(Av[k] | -Ca[k], Ca[k]) + log( 2*Ca[k] ) + log_inv_logit( A_log[k] ) + log1m_inv_logit( A_log[k] );
    target += uniform_lpdf(Bv[k] | -Cb[k], Cb[k]) + log( 2*Cb[k] ) + log_inv_logit( B_log[k] ) + log1m_inv_logit( B_log[k] );
  }
  // likelihood
  for(t in 1:T){
    rts[t,] ~ multi_normal_cholesky(mu[t,], L_H[t,]);
  }
}
//
generated quantities {
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corC;
  corr_matrix[nt] corH[T];
// Params for prediction
  vector[nt] rts_p[ahead];
  vector[nt] mu_p[ahead];
  matrix[nt,nt] rr_p[ahead];
  cov_matrix[nt] H_p[ahead];
  cholesky_factor_cov[nt] L_H_p[ahead];
  vector[2] rev_p = [0,0]';
//
//Const = multiply_lower_tri_self_transpose(Cnst);
  corC = cov2cor(Cnst);
// retrodict
  for (t in 1:T) {
    rts_out[,t] = multi_normal_rng(mu[t,], H[t,]);
       corH[t,] = cov2cor(H[t,]);
     log_lik[t] = multi_normal_lpdf(rts[t,] | mu[t,], H[t,]);
  }
// Forecast
   mu_p[1,] =  phi0 + phi * rts[T, ] + (rts[T, ]-mu[T,]) - theta * (rts[T, ]-mu[T,]);
   rr_p[1,] = ( rts[T,] - mu[T,] )*transpose( rts[T,] - mu[T,] );
    H_p[1,] = Cnst + transpose(A)*rr_p[1,]*A + transpose(B)*H[T,]*B ;    
  L_H_p[1,] = cholesky_decompose(H_p[1,]);
  rts_p[1,] = multi_normal_cholesky_rng(mu_p[1,], L_H_p[1,]);
    if(ahead >= 2) {
      for ( p in 2:ahead) {
        rev_p[2] = rts_p[p-1, 1];
        rev_p[1] = rts_p[p-1, 2];
	mu_p[p,] =  phi0 + phi * rts_p[p - 1, ] + ( rts_p[p - 1, ] - mu_p[p-1] ) - theta * ( rts_p[p - 1, ] - mu_p[p-1] );
        rr_p[p,] = ( rts_p[p - 1,] - mu_p[p - 1,] )*transpose( rts_p[p - 1,] - mu_p[p - 1,] );
         H_p[p,] = Cnst + transpose(A)*rr_p[p,]*A + transpose(B)*H_p[p-1,]*B ;  
       L_H_p[p,] = cholesky_decompose(H_p[p,]);
       rts_p[p,] = multi_normal_cholesky_rng(mu_p[p,], L_H_p[p,]);
     }
  }
}
