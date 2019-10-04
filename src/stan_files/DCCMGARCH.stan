// DCC-Parameterization
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

  // GARCH h parameters on variance metric
  vector<lower=0>[nt] c_h; 
  vector<lower=0 >[nt] a_h[Q];
  vector<lower=0, upper = 1 >[nt] b_h[P]; // TODO actually: 1 - a_h, across all Q and P...
  // GARCH q parameters 
  real<lower=0, upper = 1 > a_q; // do the mudulus thing here
  real<lower=0, upper = (1 - a_q) > b_q; //
  corr_matrix[nt] S;  // DCC keeps this constant as it is obtained from step 1. This will probably not work here (i.e. will be overwritten T times)
  // Qr1 init
  cov_matrix[nt] Qr1_init;
  // D1 init
  vector[nt] D1_init;
  // u1 init
  vector[nt] u1_init;

  real< lower = 2 > nu; // nu for student_t

  // predictor for diag variance in D
  vector[ xH_marker >= 1 ? nt : 0 ] beta;  
}
transformed parameters {
  cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  corr_matrix[nt] R[T];
  vector[nt] rr[T-1];
  vector[nt] mu[T];
  vector[nt] D[T];
  cov_matrix[nt] Qr[T];
  vector[nt] Qr_sdi[T];
  vector[nt] u[T];
  real<lower = 0> vd[nt];
  real<lower = 0> ma_d[nt];
  real<lower = 0> ar_d[nt];  

  // Initialize t=1
  mu[1,] = phi0;
  u[1,] = u1_init;
  D[1,] = D1_init;
  Qr[1,] = Qr1_init;
  L_H[1,] = cholesky_decompose(Qr[1,]);
  H[1] = L_H[1]*L_H[1]';
  R[1] = diag_matrix(rep_vector(1.0, nt));
  Qr_sdi[1] = rep_vector(1.0, nt);

  // iterations geq 2
  for (t in 2:T){
// Meanstructure model:
#include /model_components/mu.stan

    for(d in 1:nt){
      vd[d] = 0.0;
      ma_d[d] = 0.0;
      ar_d[d] = 0.0;
      // GARCH MA component
      for (q in 1:min( t-1, Q) ) {
	rr[t-q, d] = square( rts[t-q, d] - mu[t-q, d] );
	ma_d[d] = ma_d[d] + a_h[q, d]*rr[t-q, d] ;
      }
      // GARCH AR component
      for (p in 1:min( t-1, P) ) {
	ar_d[d] = ar_d[d] + b_h[p, d]*D[t-p, d];
      }

      // Predictor on diag (given in xH)
      if ( xH_marker >= 1) {
      vd[d] = c_h[d] + beta[d] * xH[t, d] + ma_d[d] + ar_d[d];
      } else if ( xH_marker == 0) {
      	vd[d] = c_h[d]  + ma_d[d] + ar_d[d];
      }

      D[t, d] = sqrt( vd[d] );
    }
    u[t,] = diag_matrix(D[t,]) \ (rts[t,]- mu[t,]); // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a
    Qr[t,] = (1 - a_q - b_q) * S + a_q * (u[t-1,] * u[t-1,]') + b_q * Qr[t-1,]; // S and UU' define dimension of Qr
    Qr_sdi[t,] = 1 ./ sqrt(diagonal(Qr[t,])); // inverse of diagonal matrix of sd's of Qr
    //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
    R[t,] = quad_form_diag(Qr[t,], Qr_sdi[t,]); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
    H[t,] = quad_form_diag(R[t,],     D[t,]);  // H = DRD; 
    L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  to_vector(beta) ~ normal(0, 3);
  // Prior for initial state
  Qr1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  to_vector(D1_init) ~ lognormal(0, 1);
  to_vector(u1_init) ~ lognormal(0, 1);
  // Prior on nu for student_t
  if ( distribution == 1 )
    nu ~ normal( nt, 50 );
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  //  to_vector(a_h) ~ normal(0, .5);
  //to_vector(b_h) ~ normal(0, .5);
  S ~ lkj_corr(nt);

  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      rts[t,] ~ multi_normal_cholesky(mu[t,], L_H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      rts[t,] ~ multi_student_t(nu, mu[t,], L_H[t,]*L_H[t,]');
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
/*   cholesky_factor_cov[nt] L_H_p[ahead]; */
/*   vector[2] rev_p = [0,0]'; */
/*   cov_matrix[nt] Qr_p[ahead]; */
/*   vector[nt] Qr_sdi_p[ahead]; */
/*   corr_matrix[nt] R_p[ahead]; */
/*   vector[nt] u_p[ahead-1]; */

  // retrodict
#include /generated/retrodict_H.stan

  /* // Forecast */
  /* mu_p[1,] =  phi0 + phi * rts[T, ] +  theta * (rts[T, ]-mu[T,]); */
  /* for(d in 1:nt){ */
  /*     rr_p[1, d] = square( rts[T, d] - mu[T, d] ); */
  /*      D_p[1, d] = sqrt( c_h[d] + a_h[d]*rr_p[1, d] +  b_h[d]*D[T,d] ); */
  /*   } */
  /* //  u_p[1,] = diag_matrix(D_p[1,]) \ (rts[t,]- mu[t,]); */
  /* Qr_p[1,] = (1 - a_q - b_q) * S + a_q * (u[T,] * u[T,]') + b_q * Qr[T,]; */
  /* Qr_sdi_p[1,] = 1 ./ sqrt(diagonal(Qr_p[1,])); */
  /* R_p[1,] = quad_form_diag(Qr_p[1,], Qr_sdi_p[1,]); */
  /* H_p[1,] = quad_form_diag(R_p[1], D_p[1]); // diag(D)*R*diag(D) */
  /* L_H_p[1,] = cholesky_decompose(H_p[1,]); */
  /* if ( distribution == 0 ) { */
  /*   rts_p[1,] = multi_normal_cholesky_rng(mu_p[1,], L_H_p[1,]); */
  /* } else if ( distribution == 1 ) { */
  /*   rts_p[1,] = multi_student_t_rng(nu, mu_p[1,], L_H_p[1,]*L_H_p[1,]'); */
  /* } */
  /* // */
  /* if(ahead >= 2) { */
  /*   for ( p in 2:ahead) { */
  /*     rev_p[2] = rts_p[p-1, 1]; */
  /*     rev_p[1] = rts_p[p-1, 2]; */
  /*     mu_p[p,] =  phi0 + phi * rts_p[p - 1, ] +  theta * ( rts_p[p - 1, ] - mu_p[p-1] ); */
  /*     for(d in 1:nt){ */
  /* 	rr_p[p, d] = square( rts_p[p-1, d] - mu_p[p-1, d] ); */
  /* 	D_p[p, d] = sqrt( c_h[d] + a_h[d]*rr_p[p-1, d] +  b_h[d]*D_p[p-1,d] ); */
  /*     } */
  /*     u_p[p-1,] = diag_matrix(D_p[p-1,]) \ (rts_p[p-1,]- mu_p[p-1,]); */
  /*     Qr_p[p,] = (1 - a_q - b_q) * S + a_q * (u_p[p-1,] * u_p[p-1,]') + b_q * Qr_p[p-1,]; */
  /*     Qr_sdi_p[p,] = 1 ./ sqrt(diagonal(Qr_p[p,])); */
  /*     R_p[p,] = quad_form_diag(Qr_p[p,], Qr_sdi_p[p,]); */
  /*     H_p[p,] = quad_form_diag(R_p[p], D_p[p]);  */
  /*     L_H_p[p,] = cholesky_decompose(H_p[p,]); */
  /*     if ( distribution == 0 ) { */
  /* 	rts_p[p,] = multi_normal_cholesky_rng(mu_p[p,], L_H_p[p,]); */
  /*     } else if ( distribution == 1 ) { */
  /* 	rts_p[p,] = multi_student_t_rng(nu, mu_p[p,], L_H_p[p,]*L_H_p[p,]'); */
  /*     } */
  /*   } */
  /* } */
}
