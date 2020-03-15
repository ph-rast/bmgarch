// pdBEKK-Parameterization
functions { 
#include /functions/cov2cor.stan
}

data {
#include /data/data.stan
}

transformed data {
  // Obtain mean and sd ove TS for prior in arma process phi0
  vector[nt] rts_m;
  vector[nt] rts_sd;
  // off diagonal elements
  int<lower = 1> od = ( nt*nt - nt ) / 2;
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

  vector<lower = 0, upper = 1>[nt] A_diag[Q];
  vector<lower = 0, upper = 1>[nt] B_diag[P];

  vector[ od ] A_lower[Q];
  vector[ od ] B_lower[P];
  vector[ od ] A_upper[Q];
  vector[ od ] B_upper[P];

    // H1 init
  cov_matrix[nt] H1_init; 
  real< lower = 2 > nu; // nu for student_t

}
transformed parameters {
  cov_matrix[nt] H[T];
  matrix[nt,nt] rr[T-1];
  vector[nt] mu[T];

  matrix[nt, nt] A_part = diag_matrix( rep_vector(0.0, nt));
  matrix[nt, nt] B_part = diag_matrix( rep_vector(0.0, nt));
  
  matrix[nt+1, nt] beta = append_row( beta0, diag_matrix(beta1) );
  row_vector[nt] C_sd;
  cov_matrix[nt] Cnst; // Const is symmetric, A, B, are not  

  // Construct square matrices with positive diagonals
  matrix[nt, nt] A_raw[Q]; 
  matrix[nt, nt] B_raw[P];

   for(q in 1:Q) {
    int L = 0;
    int U = 0;
    for( i in 1:nt ){
      for( j in 1:nt ){
	if ( i < j ) {
	  U = U + 1;
	  A_raw[q, i, j] = A_upper[q, U];
	} else if ( i > j ) {
	  L = L + 1;
	  A_raw[q, i, j] = A_lower[q, L];
	} else if (i == j ){
	  A_raw[q, i, j] = A_diag[q, i];
	}
      }
    }
  }
  

  for(p in 1:P) {
    int L = 0;
    int U = 0;
    for( i in 1:nt ){
      for( j in 1:nt ){
	if ( i < j ) {
	  U = U + 1;
	  B_raw[p, i, j] = B_upper[p, U];
	} else if ( i > j ) {
	  L = L + 1;
	  B_raw[p, i, j] = B_lower[p, L];
	} else if (i == j ){
	  B_raw[p, i, j] = B_diag[p, i];
	}
      }
    }
  }  
    
  // Initialize model parameters
  mu[1,] = phi0;
  H[1,] = H1_init;

  for (t in 2:T){    
    // Meanstructure model:
#include /model_components/mu.stan

    // reset A_part and B_part to zero for each iteration t
    A_part = diag_matrix( rep_vector(0.0, nt));
    B_part = diag_matrix( rep_vector(0.0, nt));
        
    for (q in 1:min( t-1, Q) ) {
      rr[t-q,] = ( rts[t-q,] - mu[t-q,] )*( rts[t-q,] - mu[t-q,] )';
      A_part = A_part + A_raw[q]' * rr[t-q,] * A_raw[q];
    }
    for (p in 1:min( t-1, P) ) {
      B_part = B_part + B_raw[p]' * H[t-p,] * B_raw[p];
    }
    if( xC_marker == 0 ) {
      C_sd = exp( beta0 ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst + A_part +  B_part;
    } else if( xC_marker >= 1) {
      C_sd = exp( append_col( 1.0, xC[t]' ) * beta ); 
      Cnst =  quad_form_diag(C_R, C_sd );
      H[t,] = Cnst  + A_part +  B_part;
    }
  }
}
model {
  // priors

  // Prior on nu for student_t
  nu ~ normal( nt, 50 );

  // Prior for initial state
  H1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );

  to_vector(theta) ~ std_normal();
  to_vector(phi) ~ std_normal(); 
  phi0 ~ multi_normal(rts_m, diag_matrix( rts_sd ) );

  to_vector(beta0) ~ std_normal();
  to_vector(beta1) ~ std_normal();
  C_R ~ lkj_corr( 1 );

  for(q in 1:Q) {
    to_vector(A_upper[q]) ~ std_normal();
    to_vector(A_lower[q]) ~ std_normal();
    to_vector(A_diag[q]) ~ uniform( 0, 1 );
  }
	  

  for(p in 1:P) {
    to_vector(B_upper[p]) ~ std_normal();
    to_vector(B_lower[p]) ~ std_normal();
    to_vector(B_diag[p]) ~ uniform( 0, 1 );
  }
	  
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
  matrix[nt, nt] A[Q] = A_raw;
  matrix[nt, nt] B[P] = B_raw;
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
