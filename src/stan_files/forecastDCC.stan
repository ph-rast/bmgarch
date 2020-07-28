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
#include /parameters/predH.stan
  // DF constant nu for student t
  real< lower = 2 > nu; 
  
   // GARCH h parameters on variance metric
  vector[nt] c_h;
  vector<lower=0 >[nt] a_h[Q];
  vector<lower=0, upper = 1 >[nt] b_h[P]; // TODO actually: 1 - a_h, across all Q and P...

  // GARCH q parameters
  real<lower=0, upper = 1 > a_q; //
  real<lower=0, upper = (1 - a_q) > b_q; //
  corr_matrix[nt] S;  // DCC keeps this constant

  // inits
  cov_matrix[nt] H[T];
  corr_matrix[nt] R[T];
  vector[nt] rr[T-1];
  vector[nt] mu[T];
  vector[nt] D[T];
  cov_matrix[nt] Qr[T];
  vector[nt] Qr_sdi[T];
  vector[nt] u[T];
}

generated quantities {
  // Define matrix for rts prediction
  vector[nt] rts_p[ahead + max(Q,P)];
  vector[nt] rts_forecasted[ahead];
  cov_matrix[nt] H_p[ahead + max(Q,P)];
  cov_matrix[nt] H_forecasted[ahead];
  corr_matrix[nt] R_p[ahead + max(Q,P)]; // 
  corr_matrix[nt] R_forecasted[ahead]; // 
  vector[nt] rr_p[ahead + max(Q,P)];
  vector[nt] mu_p[ahead + max(Q,P)];
  vector[nt] mu_forecasted[ahead];
  vector[nt] D_p[ahead + max(Q,P)];
  cov_matrix[nt] Qr_p[ahead + max(Q,P)];
  vector[nt] u_p[ahead + max(Q,P)];
  vector[nt] Qr_sdi_p[ahead + max(Q,P)];
  // log lik for LFO-CV
// only compute log_lik if it is actually requested 
  real log_lik[compute_log_lik ==1 ? ahead:0];
  
  // Placeholders
  real<lower = 0> vd_p[nt];
  real<lower = 0> ma_d_p[nt];
  real<lower = 0> ar_d_p[nt];
  

  // Populate with non-NA values to avoid Error in stan
  rts_p[1:(ahead + max(Q,P)), ] = rts[ 1:(ahead + max(Q,P)), ];
  H_p[  1:(ahead + max(Q,P)), ] = H[  1:(ahead + max(Q,P)), ];
  mu_p[ 1:(ahead + max(Q,P)), ] = mu[ 1:(ahead + max(Q,P)), ];
  rr_p[ 1:(ahead + max(Q,P)), ] = rr[ 1:(ahead + max(Q,P)), ];
  D_p[  1:(ahead + max(Q,P)), ] = D[  1:(ahead + max(Q,P)), ];
  u_p[  1:(ahead + max(Q,P)), ] = u[  1:(ahead + max(Q,P)), ];
  Qr_p[ 1:(ahead + max(Q,P)), ] = Qr[ 1:(ahead + max(Q,P)), ];
  Qr_sdi_p[ 1:(ahead + max(Q,P)), ] = Qr_sdi[ 1:(ahead + max(Q,P)), ];
  
  R_p[ 1:(ahead + max(Q,P)), ] = R[ 1:(ahead + max(Q,P)), ];

  // Obtain needed elements (depends on lags in Q and P) from mu and fill into mu_p
  rts_p[1:max(Q, P), ] = rts[(T-(max(Q,P)-1) ):T, ];
  H_p[  1:max(Q, P), ] =  H[ (T-(max(Q,P)-1) ):T, ];
  mu_p[ 1:max(Q, P), ] = mu[ (T-(max(Q,P)-1) ):T, ];
  // rr is of length T-1
  rr_p[ 1:max(Q, P), ] = rr[ (T-1-(max(Q,P)-1) ):(T-1), ];
  D_p[  1:max(Q, P), ] =  D[ (T - (max(Q,P)-1) ):T, ];
  u_p[  1:max(Q, P), ] =  u[ (T - (max(Q,P)-1) ):T, ];
  Qr_p[ 1:max(Q, P), ] = Qr[ (T - (max(Q,P)-1) ):T, ];
  Qr_sdi_p[ 1:max(Q, P), ] = Qr_sdi[ (T - (max(Q,P)-1) ):T, ];
  R_p[  1:max(Q, P), ] =  R[ (T - (max(Q,P)-1) ):T, ];
  
  // Forecast
  for (t in (max(Q, P) + 1 ):( max(Q, P) + ahead ) ){
    
    if( meanstructure == 0 ){
      mu_p[t, ] = phi0;
    } else if( meanstructure == 1 ){
      mu_p[t, ] = phi0 + phi * rts_p[t-1, ] + theta * (rts_p[t-1, ] - mu_p[t-1,]) ;
    }

    for(d in 1:nt){
      vd_p[d]   = 0.0;
      ma_d_p[d] = 0.0;
      ar_d_p[d] = 0.0;
      // GARCH MA component
      for (q in 1:min( t-1, Q) ) {
	rr_p[t-q, d] = square( rts_p[t-q, d] - mu_p[t-q, d] );
	ma_d_p[d] = ma_d_p[d] + a_h[q, d] * rr_p[t-q, d] ;
      }
      // GARCH AR component
      for (p in 1:min( t-1, P) ) {
      	ar_d_p[d] = ar_d_p[d] + b_h[p, d] * D_p[t-p, d]^2;
      }

      // Predictor on diag (given in xC)
      if ( xC_marker >= 1) {
	vd_p[d] = exp( c_h[d] + beta[d] * xC_c[t, d] ) + ma_d_p[d] + ar_d_p[d];
      } else if ( xC_marker == 0) {
      	vd_p[d] = exp( c_h[d] )  + ma_d_p[d] + ar_d_p[d];
      }

      D_p[t, d] = sqrt( vd_p[d] );
    }
    u_p[t - 1,] = diag_matrix(D_p[t - 1,]) \ (rts_p[t - 1,]- mu_p[t - 1,]);
    Qr_p[t, ] = (1 - a_q - b_q) * S + a_q * (u_p[t-1,] * u_p[t-1,]') + b_q * Qr_p[t-1,];
    Qr_sdi_p[t,] = 1 ./ sqrt(diagonal(Qr_p[t,]));
    R_p[t,] = quad_form_diag(Qr_p[t,], Qr_sdi_p[t,]);
    
    //    H_p[t,] = quad_form_diag(R_p[t,],     D_p[t,]);
    H_p[t,] = diag_matrix(D_p[t,]) * R_p[t,] * diag_matrix(D_p[t,]);
       
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
    /* 	    } */
    /* 	  } */
    /* } */
  }
  rts_forecasted = rts_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  H_forecasted = H_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  R_forecasted = R_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
  mu_forecasted = mu_p[max(Q, P) + 1 : (max(Q, P) + ahead)];
#include /generated/forecast_log_lik.stan
}
