// CCC-Parameterization
functions { 
#include /functions/cov2cor.stan
}
data {
#include /data/data.stan
}
transformed data {
  // Obtain mean and sd over TS for prior in arma process phi0                                                                                                                                                 
  vector[nt] rts_m;
  vector[nt] rts_sd;

#include /transformed_data/xh_marker.stan
   
  for ( i in 1:nt ){
    rts_m[i] = mean(rts[,i]);
    rts_sd[i] = sd(rts[,i]);
  }
}
parameters {
 // ARMA parameters 
#include /parameters/arma.stan
  // predictor for H
#include /parameters/predH.stan

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

}
transformed parameters {
  cov_matrix[nt] H[T];
  vector[nt] rr[T-1];
  vector[nt] mu[T];
  vector[nt] D[T];
  real<lower = 0> vd[nt];
  real<lower = 0> ma_d[nt];
  real<lower = 0> ar_d[nt];
  // Initialize t=1
  // Check "Order Sensitivity and Repeated Variables" in stan reference manual
  mu[1,] = phi0;
  //u[1,] = diagonal(sigma1);  
  D[1,] = D1_init;
  H[1,] = quad_form_diag(R, D[1,]);  // H = DRD; 
  // iterations geq 2
  for (t in 2:T){

    // Meanstructure model:
#include /model_components/mu.stan
    
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
	ar_d[d] = ar_d[d] + b_h[p, d]*D[t-p, d]^2;
      }
      if ( xC_marker >= 1) {
	vd[d] = exp( c_h[d] + beta[d] * xC[t, d] ) + ma_d[d] + ar_d[d];
      } else if ( xC_marker == 0) {
      	vd[d] = exp( c_h[d] )  + ma_d[d] + ar_d[d];
      }
      D[t, d] = sqrt( vd[d] );
    }
  H[t,] = quad_form_diag(R, D[t,]);  // H = DRD;
  }
}

model {
  // priors
  to_vector(beta) ~ normal(0, 1);
  to_vector(c_h) ~ normal(-2, 4);
  if ( distribution == 1 )
    nu ~ normal( nt, 50 );
  to_vector(D1_init) ~ lognormal(0, 1);
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  phi0 ~ multi_normal(rts_m, diag_matrix( rts_sd ) );
  R ~ lkj_corr( 1 );
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

  // retrodict
#include /generated/retrodict_H.stan
  
}
