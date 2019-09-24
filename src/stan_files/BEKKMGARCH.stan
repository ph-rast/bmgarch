// BEKK-Parameterization
functions { 
#include /functions/cov2cor.stan
}
data {
  int<lower=2> T;
  int<lower=1> nt;    // number of time series
  int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
  int<lower=1> P; // AR component in MGARCH(P,Q), matrix B  
  vector[nt] rts[T];  // multivariate time-series
  vector[nt] xH[T];  // time-varying predictor for conditional H
  int<lower=0> ahead; // how many ahead predictions
  int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
}
transformed data {
  matrix[nt, nt] xH_m[T];
  real<lower = 0> xH_marker = 0.0;
  real<lower = 0> cp;
  for( t in 1:T ){
    xH_m[t] = diag_matrix( xH[t] );
  // add a variable that notes if xH is null or actually a predictor
    cp = sum( xH[t] );
    if( cp != 0)
      xH_marker = xH_marker + 1;
  }
}
parameters { 
  //  cholesky_factor_cov[nt] Cnst; // Const is symmetric, A, B, are not
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  
  // construct A, so that one element (a11) can be constrained to be non-negative
  real<lower = 0, upper = 1> Ap11[Q];
  row_vector[nt-2] Ap1k[Q];
  matrix[nt-1, nt-1] Ap_sub[Q];
  //
  real<lower = 0, upper = 1> Bp11[P];
  row_vector[nt-2] Bp1k[P];
  matrix[nt-1, nt-1] Bp_sub[P];
  //
  vector[nt] A_log[Q];
  vector[nt] B_log[P];
  // ARMA parameters
  vector[nt] phi0; 
  matrix[nt,nt] phi;
  matrix[nt,nt] theta;
  // H1 init
  cov_matrix[nt] H1_init;
  real< lower = 2 > nu; // nu for student_t
  // predictor for H
  cov_matrix[nt] beta;
}
transformed parameters {
  cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];
  matrix[nt, nt] A[Q];
  matrix[nt, nt] B[P];
  vector[nt] Ca[Q]; // Upper (and lower) boundary for A 
  vector[nt] Av[Q];
  vector[nt] Cb[P]; // Upper (and lower) boundary for B
  vector[nt] Bv[P]; 
  matrix[nt, nt -1 ] Ap[Q];
  matrix[nt, nt -1 ] Bp[P];
  matrix[nt, nt] A_part = diag_matrix( rep_vector(0.0, nt));
  matrix[nt, nt] B_part = diag_matrix( rep_vector(0.0, nt));  
  for ( q in 1:Q )
    Ap[q] = append_row(append_col(Ap11[q], Ap1k[q]), Ap_sub[q]);
  for ( p in 1:P )
    Bp[p] = append_row(append_col(Bp11[p], Bp1k[p]), Bp_sub[p]);
  for ( i in 1:nt) {
    for ( q in 1:Q ){
      Ca[q,i] = sqrt( 1 - dot_self(Ap[q,i]) );
      Av[q,i] = -Ca[q,i] + 2*Ca[q,i] * inv_logit( A_log[q,i] );
    }
    for ( p in 1:P ){ 
      Cb[p,i] = sqrt( 1 - dot_self(Bp[p,i]) );
      Bv[p,i] = -Cb[p,i] + 2*Cb[p,i] * inv_logit( B_log[p,i] );
    }
  }
  
  for ( q in 1:Q ) 
    A[q] = append_col(Ap[q], Av[q]);

  for ( p in 1:P )
    B[p] = append_col(Bp[p], Bv[p]);

  // Initialize
  H[1,] = H1_init;
  L_H[1,] = cholesky_decompose(H[1,]); // cf. p 69 in stan manual for how to index
  // Means AR + MA 
  mu[1,] = phi0 + phi * rts[1, ] + theta * (rts[1, ] - phi0) ;

  //
  for (t in 2:T){
    mu[t, ] = phi0 + phi * rts[t-1, ] +  theta * (rts[t-1, ] - mu[t-1,]) ;
    for (q in 1:min( t-1, Q) ) {
      rr[t-q,] = ( rts[t-q,] - mu[t-q,] )*( rts[t-q,] - mu[t-q,] )';
      A_part = A_part + A[q]' * rr[t-q,] * A[q];
    }
    for (p in 1:min( t-1, P) ) {
      B_part = B_part + B[p]' * H[t-p,] * B[p];
    }
    if( xH_marker == 0 ) {
      H[t,] = Cnst + A_part +  B_part;} else {
      H[t,] = Cnst + beta * xH_m[t]  + A_part +  B_part;
    }
    L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  // Prior on nu for student_t
  if ( distribution == 1 )
    nu ~ normal( nt, 50 );
  // Prior for initial state
  H1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  Cnst ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
 
  for( k in 1:nt){ 
     for( q in 1:Q ) {
       target += uniform_lpdf(Av[q,k] | -Ca[q,k], Ca[q,k]) + log( 2*Ca[q,k] ) + log_inv_logit( A_log[q,k] ) + log1m_inv_logit( A_log[q,k] );
     }
     for ( p in 1:P ) {
       target += uniform_lpdf(Bv[p,k] | -Cb[p,k], Cb[p,k]) + log( 2*Cb[p,k] ) + log_inv_logit( B_log[p,k] ) + log1m_inv_logit( B_log[p,k] );
     }
  }
  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      //      rts[t,] ~ multi_normal_cholesky(mu[t,], L_H[t,]);
      target += multi_normal_cholesky_lpdf( rts[t, ] | mu[t,], L_H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      // rts[t,] ~ multi_student_t(nu, mu[t,], L_H[t,]*L_H[t,]');
      target += multi_student_t_lpdf( rts[t, ] | nu, mu[t,], L_H[t,]*L_H[t,]');
    }
  }
}
//
generated quantities {
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corC;
  corr_matrix[nt] corH[T];
/* // Params for prediction */
/*   vector[nt] rts_p[ahead]; */
/*   vector[nt] mu_p[ahead]; */
/*   matrix[nt,nt] rr_p[ahead]; */
/*   cov_matrix[nt] H_p[ahead]; */
/*   cholesky_factor_cov[nt] L_H_p[ahead]; */
/*   vector[2] rev_p = [0,0]'; */
//
//Const = multiply_lower_tri_self_transpose(Cnst);
  corC = cov2cor(Cnst);

  // retrodict
#include /generated/retrodict_H.stan

  /* // Forecast */
  /*  mu_p[1,] =  phi0 + phi * rts[T, ] +  theta * (rts[T, ]-mu[T,]); */
  /*  rr_p[1,] = ( rts[T,] - mu[T,] )*transpose( rts[T,] - mu[T,] ); */
  /*   H_p[1,] = Cnst + transpose(A)*rr_p[1,]*A + transpose(B)*H[T,]*B ;     */
  /* L_H_p[1,] = cholesky_decompose(H_p[1,]); */
  /* if ( distribution == 0 ) { */
  /*   rts_p[1,] = multi_normal_cholesky_rng(mu_p[1,], L_H_p[1,]); */
  /* } else if ( distribution == 1 ) { */
  /*   rts_p[1,] = multi_student_t_rng(nu, mu_p[1,], L_H_p[1,]*L_H_p[1,]'); */
  /* } */
  /* // */
  /*   if(ahead >= 2) { */
  /*     for ( p in 2:ahead) { */
  /*       rev_p[2] = rts_p[p-1, 1]; */
  /*       rev_p[1] = rts_p[p-1, 2]; */
  /* 	mu_p[p,] =  phi0 + phi * rts_p[p - 1, ] +  theta * ( rts_p[p - 1, ] - mu_p[p-1] ); */
  /*       rr_p[p,] = ( rts_p[p - 1,] - mu_p[p - 1,] )*transpose( rts_p[p - 1,] - mu_p[p - 1,] ); */
  /*        H_p[p,] = Cnst + transpose(A)*rr_p[p,]*A + transpose(B)*H_p[p-1,]*B ;   */
  /*      L_H_p[p,] = cholesky_decompose(H_p[p,]); */
  /*      if ( distribution == 0 ) { */
  /* 	rts_p[p,] = multi_normal_cholesky_rng(mu_p[p,], L_H_p[p,]); */
  /*     } else if ( distribution == 1 ) { */
  /* 	rts_p[p,] = multi_student_t_rng(nu, mu_p[p,], L_H_p[p,]*L_H_p[p,]'); */
  /*     } */
  /*    } */
  /* } */
}
