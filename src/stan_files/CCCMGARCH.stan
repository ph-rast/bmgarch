// CCC-Parameterization
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
  int<lower=1> nt;                    // number of time series
  vector[nt] rts[T];  // multivariate time-series
  int<lower=0> ahead; // how many ahead predictions 
  cov_matrix[nt] sigma1;
}
transformed data {
  // Reverse the rts vector
  vector[nt] rev[T];
  for( i in 1:nt ) {
    rev[,i] = rts[,nt-i+1];
  }  
}
parameters { 
  vector[nt] phi0; 
  matrix[nt,nt] phi;
  matrix[nt,nt] theta;
  // GARCH h parameters on variance metric
  vector<lower=0>[nt] c_h; 
  vector<lower=0 >[nt] a_h;
  vector<lower=0 >[nt] b_h; // not sure if this upper def works with vectors
  // GARCH q parameters 
//  real<lower=0, upper = 1 > a_q; // do the mudulus thing here
//  real<lower=0, upper = (1 - a_q) > b_q; //
  corr_matrix[nt] R;  
}
transformed parameters {
  //cholesky_factor_cov[nt] L_H[T];
  cov_matrix[nt] H[T];
  //corr_matrix[nt] R;
  vector[nt] rr[T-1];
  vector[nt] mu[T];
  vector[nt] D[T];
//  cov_matrix[nt] Q[T];
//  vector[nt] Q_sdi[T];
//  vector[nt] u[T];
  // Initialize t=1
    mu[1,] = phi0 + phi * rts[1, ] + (rts[1, ] - phi0) - theta * (rts[1, ] - phi0) ;
  
  //u[1,] = diagonal(sigma1);
  D[1,] = diagonal(sigma1);
  H[1,] = quad_form_diag(R,     D[1,]);  // H = DRD; 
  //Q[1,] = sigma1;
//  L_H[1,] = cholesky_decompose(Q[1,]);
//  H[1] = L_H[1]*L_H[1]';
//  R[1] = diag_matrix(rep_vector(1.0, nt));
//  Q_sdi[1] = rep_vector(1.0, nt);
  // iterations geq 2
  for (t in 2:T){
    mu[t, ] = phi0 + phi * rts[t-1, ] + (rts[t-1, ] - mu[t-1,]) - theta * (rts[t-1, ] - mu[t-1,]) ;
    for(d in 1:nt){ 
      rr[t-1,d] = square( rts[t-1,d] - mu[t-1,d] );
      D[t,d] = sqrt( c_h[d] + a_h[d]*rr[t-1,d] +  b_h[d]*D[t-1,d] );  
    }
//    u[t,] = diag_matrix(D[t,]) \ (rts[t,]- mu[t,]); // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a
//    Q[t,] = (1 - a_q - b_q) * S + a_q * (u[t-1,] * u[t-1,]') + b_q * Q[t-1,]; // S and UU' define dimension of Q
//    Q_sdi[t,] = 1 ./ sqrt(diagonal(Q[t,])); // inverse of diagonal matrix of sd's of Q
    //    R[t,] = quad_form_diag(Q[t,], inv(sqrt(diagonal(Q[t,]))) ); // Q_sdi[t,] * Q[t,] * Q_sdi[t,];
//    R = quad_form_diag(Q[t,], Q_sdi[t,]); // Q_sdi[t,] * Q[t,] * Q_sdi[t,];
    H[t,] = quad_form_diag(R,     D[t,]);  // H = DRD; 
//    L_H[t,] = cholesky_decompose(H[t,]);
  }
}
model {
  // priors
  to_vector(theta) ~ normal(0, 1);
  to_vector(phi) ~ normal(0, 1);
  to_vector(phi0) ~ normal(0, 1);
  to_vector(a_h) ~ normal(0, .5);
  to_vector(b_h) ~ normal(0, .5);
  R ~ lkj_corr(nt);
  // likelihood
  for(t in 1:T){
    rts[t,] ~ multi_normal(mu[t,], H[t,]);
  }
}
generated quantities {
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corH[T];
// Params for prediction
  vector[nt] rts_p[ahead];
  vector[nt] mu_p[ahead];
  vector[nt] rr_p[ahead];
  vector[nt] D_p[ahead];
  cov_matrix[nt] H_p[ahead];
 // cholesky_factor_cov[nt] L_H_p[ahead];
  vector[2] rev_p = [0,0]';
//  cov_matrix[nt] Q_p[ahead];
//  vector[nt] Q_sdi_p[ahead];
//  corr_matrix[nt] R_p[ahead];
//  vector[nt] u_p[ahead-1];
// retrodict
  for (t in 1:T) {
    rts_out[,t] = multi_normal_rng(mu[t,], H[t,]);
       corH[t,] = cov2cor(H[t,]);
     log_lik[t] = multi_normal_lpdf(rts[t,] | mu[t,], H[t,]);
  }
// Forecast
  mu_p[1,] =  phi0 + phi * rts[T, ] + (rts[T, ]-mu[T,]) - theta * (rts[T, ]-mu[T,]);
  for(d in 1:nt){
      rr_p[1, d] = square( rts[T, d] - mu[T, d] );
       D_p[1, d] = sqrt( c_h[d] + a_h[d]*rr_p[1, d] +  b_h[d]*D[T,d] );
    }
  //  u_p[1,] = diag_matrix(D_p[1,]) \ (rts[t,]- mu[t,]);
//  Q_p[1,] = (1 - a_q - b_q) * S + a_q * (u[T,] * u[T,]') + b_q * Q[T,];
//  Q_sdi_p[1,] = 1 ./ sqrt(diagonal(Q_p[1,]));
//  R_p[1,] = quad_form_diag(Q_p[1,], Q_sdi_p[1,]);
  H_p[1,] = quad_form_diag(R, D_p[1]); // diag(D)*R*diag(D)
  //L_H_p[1,] = cholesky_decompose(H_p[1,]);
  rts_p[1,] = multi_normal_rng(mu_p[1,], H_p[1,]);
  //
  if(ahead >= 2) {
    for ( p in 2:ahead) {
      rev_p[2] = rts_p[p-1, 1];
      rev_p[1] = rts_p[p-1, 2];
      mu_p[p,] =  phi0 + phi * rts_p[p - 1, ] + ( rts_p[p - 1, ] - mu_p[p-1] ) - theta * ( rts_p[p - 1, ] - mu_p[p-1] );
      for(d in 1:nt){
	rr_p[p, d] = square( rts_p[p-1, d] - mu_p[p-1, d] );
	D_p[p, d] = sqrt( c_h[d] + a_h[d]*rr_p[p-1, d] +  b_h[d]*D_p[p-1,d] );
      }
//      u_p[p-1,] = diag_matrix(D_p[p-1,]) \ (rts_p[p-1,]- mu_p[p-1,]);
//      Q_p[p,] = (1 - a_q - b_q) * S + a_q * (u_p[p-1,] * u_p[p-1,]') + b_q * Q_p[p-1,];
//      Q_sdi_p[p,] = 1 ./ sqrt(diagonal(Q_p[p,]));
//      R_p[p,] = quad_form_diag(Q_p[p,], Q_sdi_p[p,]);
      H_p[p,] = quad_form_diag(R, D_p[p]); 
    //  L_H_p[p,] = cholesky_decompose(H_p[p,]);
      rts_p[p,] = multi_normal_rng(mu_p[p,], H_p[p,]);
    }
  }
}
