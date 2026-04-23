// Constant covariance parameterization (no ARCH/GARCH dynamics)
// Baseline model for LFO-CV comparison with MGARCH parameterizations.
// H is a fixed positive-definite covariance matrix, parameterized via constant
// log-variances (c_h) and a constant LKJ correlation matrix (R).
// Supports constant, ARMA(1,1), and VAR(1) mean structures.
functions {
#include /functions/cov2cor.stan
#include /functions/jacobian.stan
}
data {
#include /data/data.stan
}
transformed data {
  vector[nt] rts_m;
  vector[nt] rts_sd;

#include /transformed_data/xh_marker.stan

  if( meanstructure == 0 ){
    for ( i in 1:nt ){
      rts_m[i] = mean(rts[,i]);
      rts_sd[i] = sd(rts[,i]);
    }
  } else if (meanstructure == 1 || meanstructure == 2){
    for ( i in 1:nt ){
      rts_m[i] = rts[1,i];
      rts_sd[i] = sd(rts[,i]);
    }
  }
}
parameters {
#include /parameters/arma.stan
#include /parameters/predH.stan

  // Constant log-variances (variance = exp(c_h))
  vector[nt] c_h;
  // Constant correlation matrix
  corr_matrix[nt] R;
  // Degrees of freedom for Student-t (used when distribution == 1)
  real<lower=2> nu;
}
transformed parameters {
  cov_matrix[nt] H_const;
  array[T] vector[nt] mu;

  H_const = quad_form_diag(R, sqrt(exp(c_h)));

  mu[1,] = phi0;
  for (t in 2:T) {
#include /model_components/mu.stan
  }
}
model {
  to_vector(c_h)   ~ std_normal();
  to_vector(beta)  ~ std_normal();
  to_vector(theta) ~ std_normal();
  to_vector(phi)   ~ std_normal();
  phi0 ~ multi_normal(rts_m, diag_matrix(rts_sd));
  R ~ lkj_corr(1);
  if (distribution == 1)
    nu ~ normal(nt, 50);

  if (distribution == 0) {
    for (t in 1:T)
      rts[t,] ~ multi_normal(mu[t,], H_const);
  } else if (distribution == 1) {
    for (t in 1:T)
      rts[t,] ~ multi_student_t(nu, mu[t,], H_const);
  }
}
generated quantities {
  matrix[nt,T] rts_out;
  array[T] real log_lik;
  // H and corH replicated T times for API consistency with fitted.bmgarch / plot.
  // Declared as matrix (not cov_matrix/corr_matrix) to avoid rep_array
  // constraint-checking issues under parallel chain execution.
  array[T] matrix[nt,nt] H;
  array[T] matrix[nt,nt] corH;
  {
    matrix[nt,nt] corH_const = cov2cor(H_const);
    for (t in 1:T) {
      H[t]    = H_const;
      corH[t] = corH_const;
    }
  }
  vector<lower=0>[nt] c_h_var = exp(c_h);

  if (distribution == 0) {
    for (t in 1:T) {
      rts_out[,t] = multi_normal_rng(mu[t,], H_const);
      log_lik[t]  = multi_normal_lpdf(rts[t,] | mu[t,], H_const);
    }
  } else if (distribution == 1) {
    for (t in 1:T) {
      rts_out[,t] = multi_student_t_rng(nu, mu[t,], H_const);
      log_lik[t]  = multi_student_t_lpdf(rts[t,] | nu, mu[t,], H_const);
    }
  }
}
