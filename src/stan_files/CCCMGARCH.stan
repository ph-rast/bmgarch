// CCC-Parameterization
functions { 
#include /functions/cov2cor.stan
}
data {
  int<lower=2> T;
  int<lower=1> nt;                    // number of time series
  int<lower=1> Q; // MA component in MGARCH(P,Q), vector A
  int<lower=1> P; // AR component in MGARCH(P,Q), vector B
  vector[nt] rts[T];  // multivariate time-series
  vector[nt] xH[T];  // time-varying predictor for conditional H
  int<lower=0> ahead; // how many ahead predictions
  int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
  int<lower=0, upper=1> meanstructure; // Select model for location
}
transformed data {
#include /transformed_data/xh_marker.stan  
}
parameters {
  // ARMA parameters
  vector[nt] phi0; 
  matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] phi;
  matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] theta;
  // GARCH h parameters on variance metric
  vector<lower=0>[nt] c_h; 
  vector<lower=0, upper = 1 >[nt] a_h[Q];
  vector<lower=0, upper = 1 >[nt] b_h[P]; // TODO actually: 1 - a_h, across all Q and P...
  // GARCH constant correlation 
  corr_matrix[nt] R;
  // D1 init
  vector<lower = 0>[nt] D1_init;
  // DF constant nu for student t
  real< lower = 2 > nu;
  // predictor for H 
  vector[ xH_marker >= 1 ? nt : 0 ] beta; 
}
transformed parameters {
  //cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  //corr_matrix[nt] R;
  vector[nt] rr[T-1];
  vector[nt] mu[T];
  vector[nt] D[T];
  real<lower = 0> vd[nt];
  real<lower = 0> ma_d[nt];
  real<lower = 0> ar_d[nt];
  // Initialize t=1
  // Check "Order Sensitivity and Repeated Variables" in stan reference manual
  mu[1,] = phi0; // + phi * rts[1, ] + theta * (rts[1, ] - phi0) ;
  //u[1,] = diagonal(sigma1);  
  D[1,] = D1_init;
  H[1,] = quad_form_diag(R, D[1,]);  // H = DRD; 
  // iterations geq 2
  for (t in 2:T){

    // location:
    if( meanstructure == 1 ){
    mu[t, ] = phi0 + phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
    } else if ( meanstructure == 0 ){
      mu[t, ] = phi0;
    }
    
    // scale: SD's of D
    for (d in 1:nt) {
      vd[d] = 0.0;
      ma_d[d] = 0.0;
      ar_d[d] = 0.0;
      // MA component
      for (q in 1:min( t-1, Q) ) {
	rr[t-q, d] = square( rts[t-q, d] - mu[t-q, d] );
	ma_d[d] = ma_d[d] + a_h[q, d]*rr[t-q, d] ;
      }
      for (p in 1:min( t-1, P) ) {
	ar_d[d] = ar_d[d] + b_h[p, d]*D[t-p, d];
      }
      if ( xH_marker >= 1) {
      vd[d] = c_h[d] + beta[d] * xH[t, d] + ma_d[d] + ar_d[d];
      } else if ( xH_marker == 0) {
      	vd[d] = c_h[d]  + ma_d[d] + ar_d[d];
      }
      D[t, d] = sqrt( vd[d] );
    }
  H[t,] = quad_form_diag(R, D[t,]);  // H = DRD;
  }
}

model {
  // priors
  to_vector(beta) ~ normal(0, 3);
  if ( distribution == 1 )
    nu ~ normal( nt, 50 );
  to_vector(D1_init) ~ lognormal(0, 1);
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  R ~ lkj_corr(nt);
  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      //      rts[t,] ~ multi_normal(mu[t,], H[t,]);
      target += multi_normal_lpdf( rts[t, ] | mu[t,], H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      //      rts[t,] ~ multi_student_t(nu, mu[t,], H[t,]);
      target += multi_student_t_lpdf( rts[t, ] | nu, mu[t,], H[t,]);
    }
  }
}
generated quantities {
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corH[T];
/* // Params for prediction */
/*   vector[nt] rts_p[ahead]; */
/*   vector[nt] mu_p[ahead]; */
/*   vector[nt] rr_p[ahead]; */
/*   vector[nt] D_p[ahead]; */
/*   cov_matrix[nt] H_p[ahead]; */
/*   vector[2] rev_p = [0,0]'; */

// retrodict
#include /generated/retrodict_H.stan
  
/* // Forecast */
/*   mu_p[1,] =  phi0 + phi * rts[T, ] +  theta * (rts[T, ]-mu[T,]); */
/*   for(d in 1:nt){ */
/*       rr_p[1, d] = square( rts[T, d] - mu[T, d] ); */
/*        D_p[1, d] = sqrt( c_h[d] + a_h[d]*rr_p[1, d] +  b_h[d]*D[T,d] ); */
/*     } */
/*   H_p[1,] = quad_form_diag(R, D_p[1]); */
/*   if ( distribution == 0 ) { */
/*     rts_p[1,] = multi_normal_rng(mu_p[1,], H_p[1,]); */
/*   } else if ( distribution == 1 ) { */
/*     rts_p[1,] = multi_student_t_rng(nu, mu_p[1,], H_p[1,]); */
/*   } */
/*   // */
/*   if(ahead >= 2) { */
/*     for ( p in 2:ahead) { */
/*       rev_p[2] = rts_p[p-1, 1]; */
/*       rev_p[1] = rts_p[p-1, 2]; */
/*       mu_p[p,] =  phi0 + phi * rts_p[p - 1, ] + theta * ( rts_p[p - 1, ] - mu_p[p-1] ); */
/*       for(d in 1:nt){ */
/* 	rr_p[p, d] = square( rts_p[p-1, d] - mu_p[p-1, d] ); */
/* 	D_p[p, d] = sqrt( c_h[d] + a_h[d]*rr_p[p-1, d] +  b_h[d]*D_p[p-1,d] ); */
/*       } */
/*       H_p[p,] = quad_form_diag(R, D_p[p]); */
/*       if ( distribution == 0 ) { */
/* 	rts_p[p,] = multi_normal_rng(mu_p[p,], H_p[p,]); */
/*       } else if ( distribution == 1 ) { */
/* 	rts_p[p,] = multi_student_t_rng(nu, mu_p[p,], H_p[p,]); */
/*       } */
/*     } */
/*   } */
}
