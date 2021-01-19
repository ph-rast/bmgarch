// Constant variance - No garch
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

  if( meanstructure == 0 ){
    for ( i in 1:nt ){
      rts_m[i] = mean(rts[,i]);
      rts_sd[i] = sd(rts[,i]);
    }
  } else if (meanstructure == 1 ){
      // set rts_m to first element in ts
      for ( i in 1:nt ){
	rts_m[i] = rts[1,i];
	rts_sd[i] = sd(rts[,i]);
      }
    }
}

parameters { 
  // ARMA parameters 
#include /parameters/arma.stan

  // predictor for H
  //cov_matrix[ xC_marker >= 1 ? nt : 0 ] beta;
  row_vector[nt] beta0;
  vector[nt] beta1;
  
  // Separate into C_sd*C_R*C_sd
  // C_sd is defined in tp, as function of betas
  corr_matrix[nt] C_R;

    // H1 init
  cov_matrix[nt] H1_init; 
  real< lower = 2 > nu; // nu for student_t

}
transformed parameters {
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];


  matrix[nt+1, nt] beta = append_row( beta0, diag_matrix(beta1) );
  row_vector[nt] C_sd;
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  
    
  // Initialize model parameters
  mu[1,] = phi0;
  H[1,] = H1_init;

  for (t in 2:T){    
    // Meanstructure model:
#include /model_components/mu.stan

    if( xC_marker == 0 ) {
      C_sd = exp( beta0 ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst ;
    } else if( xC_marker >= 1) {
      C_sd = exp( append_col( 1.0, xC[t]' ) * beta ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst ;
    }
  }
}
model {
  // priors
  // https://mc-stan.org/documentation/case-studies/mle-params.html

  // Prior on nu for student_t
  nu ~ normal( nt, 50 );

  // Prior for initial state
  H1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );

  to_vector(theta) ~ std_normal();
  to_vector(phi) ~ std_normal(); 
  //  to_vector(phi0) ~ normal();
  phi0 ~ multi_normal(rts_m, diag_matrix( rts_sd ) );

  to_vector(beta0) ~ std_normal();
  to_vector(beta1) ~ std_normal();
  C_R ~ lkj_corr( 1 );
   
  // likelihood
  if ( distribution == 0 ) {
    for(t in 1:T){
      target += multi_normal_lpdf(rts[t, ]| mu[t,], H[t,]);
    }
  } else if ( distribution == 1 ) {
    for(t in 1:T){
      target += multi_student_t_lpdf( rts[t, ] | nu, mu[t,], H[t,]);
    }
  }
}
//
generated quantities {
  matrix[nt,T] rts_out;
  real log_lik[T];
  corr_matrix[nt] corC;
  corr_matrix[nt] corH[T];
  row_vector[nt] C_var;

//Const = multiply_lower_tri_self_transpose(Cnst);
  corC = C_R;
  C_var = exp(2*beta0 );

  // retrodict
#include /generated/retrodict_H.stan
}
