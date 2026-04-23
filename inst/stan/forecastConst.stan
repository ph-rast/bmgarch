// Forecast generated quantities for the constant covariance model.
// H_const is reconstructed from posterior draws of c_h and R; no GARCH
// propagation is needed. Mean structure is propagated forward via VAR/ARMA/const.
data {
#include /data/gq_data.stan
}
transformed data {
  matrix[ahead + max(Q,P), nt] xC_c;
#include /transformed_data/xh_marker.stan
  for (i in 1:max(Q,P)) {
    xC_c[i] = xC[T - (max(Q,P) - 1)]';
  }
  for (i in 1:ahead) {
    xC_c[i + max(Q,P)] = xC_p[i]';
  }
}
parameters {
#include /parameters/arma.stan
#include /parameters/predH.stan
  real<lower=2> nu;
  vector[nt] c_h;
  corr_matrix[nt] R;
  // mu from training needed to seed ARMA mean structure
  array[T] vector[nt] mu;
}
generated quantities {
  cov_matrix[nt] H_const = quad_form_diag(R, sqrt(exp(c_h)));

  array[ahead + max(Q,P)] cov_matrix[nt] H_p;
  array[ahead] cov_matrix[nt] H_forecasted;
  array[ahead + max(Q,P)] corr_matrix[nt] R_p = rep_array(R, ahead + max(Q,P));
  array[ahead] corr_matrix[nt] R_forecasted = rep_array(R, ahead);

  array[ahead + max(Q,P)] vector[nt] mu_p;
  array[ahead] vector[nt] mu_forecasted;
  array[ahead + max(Q,P)] vector[nt] rts_p;
  array[ahead] vector[nt] rts_forecasted;

  array[compute_log_lik == 1 ? ahead : 0] real log_lik;

  // Seed with last max(Q,P) steps from training
  mu_p[ 1:max(Q,P), ] = mu[ (T-(max(Q,P)-1)):T, ];
  rts_p[1:max(Q,P), ] = rts[(T-(max(Q,P)-1)):T, ];
  for (i in 1:max(Q,P))
    H_p[i,] = H_const;

  // Forecast ahead steps
  for (t in (max(Q,P)+1):(max(Q,P)+ahead)) {
    if (meanstructure == 0) {
      mu_p[t,] = phi0;
    } else if (meanstructure == 1) {
      mu_p[t,] = phi0 + phi * rts_p[t-1,] + theta * (rts_p[t-1,] - mu_p[t-1,]);
    } else if (meanstructure == 2) {
      mu_p[t,] = phi0 + phi * rts_p[t-1,];
    }
    H_p[t,] = H_const;

#include /generated/forecast_sampling.stan
  }

  rts_forecasted = rts_p[max(Q,P)+1:(max(Q,P)+ahead)];
  H_forecasted   = H_p[  max(Q,P)+1:(max(Q,P)+ahead)];
  mu_forecasted  = mu_p[ max(Q,P)+1:(max(Q,P)+ahead)];

#include /generated/forecast_log_lik.stan
}
