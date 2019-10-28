data {

  // Data 
  int<lower=2> T;  // Original TS length
  int<lower=0> ahead; // how many ahead predictions

  int<lower=1> nt;    // number of time series
  int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
  int<lower=1> P; // AR component in MGARCH(P,Q), matrix B

  vector[nt] rts[T];  // multivariate time-series 
  vector[nt] xH[T];  // time-varying predictor for conditional H

  int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
  int<lower=0, upper=1> meanstructure; // Select model for location

  // Parameters from fitted object
  vector[nt] phi0; 
  matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] phi;
  matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] theta;

  vector[nt] mu[T];  
  vector[nt] rr[T-1]; 
  vector[nt] D[T]; 
  cov_matrix[nt] H[T];

  real<lower = nt> nu;

}
parameters {
}
model {

}
generated quantities {
  // Params for prediction
  vector[nt] D_p[ahead];
  cov_matrix[nt] H_p[ahead];
  real<lower = 0> vd[nt];
  real<lower = 0> ma_d[nt]; 
  real<lower = 0> ar_d[nt];

  // Define Vector that contains max of Q or P lag plus the forecasted ahed
  vector[nt] mu_p[ahead + max(Q,P)];

  // Obtain needed elements from mu and fill into mu_p
  mu_p[1:max(Q, P)] = mu[T-max(Q,P):T];
  
  // Forecast
  for (t in (max(Q, P) + 1 ):( max(Q, P) + ahead ) ){

    if( meanstructure == 0 ){
      mu_p[t, ] = phi0;
    } else if( meanstructure == 1 ){
      mu_p[t, ] = phi0 + phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
    }
 
    /* // scale: SD's of D */
    /* for (d in 1:nt) { */
    /*   vd[d] = 0.0; */
    /*   ma_d[d] = 0.0; */
    /*   ar_d[d] = 0.0; */
    /*   // MA component */
    /*   for (q in 1:min( t-1, Q) ) { */
    /* 	rr[t-q, d] = square( rts[t-q, d] - mu[t-q, d] ); */
    /* 	ma_d[d] = ma_d[d] + a_h[q, d]*rr[t-q, d] ; */
    /*   } */
    /*   for (p in 1:min( t-1, P) ) { */
    /* 	ar_d[d] = ar_d[d] + b_h[p, d]*D[t-p, d]; */
    /*   } */
    /*   if ( xH_marker >= 1) { */
    /* 	vd[d] = c_h[d] + beta[d] * xH[t, d] + ma_d[d] + ar_d[d]; */
    /*   } else if ( xH_marker == 0) { */
    /*   	vd[d] = c_h[d]  + ma_d[d] + ar_d[d]; */
    /*   } */
    /*   D[t, d] = sqrt( vd[d] ); */
    /* } */
    /* H[t,] = quad_form_diag(R, D[t,]);  // H = DRD; */

    /* // likelihood */
    /* if ( distribution == 0 ) { */
    /*   rts[t,] = multi_normal_rng( mu[t,], H[t,]); */
    /* } else if ( distribution == 1 ) { */
    /*   rts[t,] = multi_student_t_rng( nu, mu[t,], H[t,]); */
    /* } */
  }
}

