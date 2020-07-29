real a_b_scale_jacobian(real a, real b, real value) {
  real invlogit_value = inv_logit(value);
  real out = log(b - a) + log(invlogit_value) + log1m(invlogit_value);
  return(out);
}

real a_b_scale(real a, real b, real value) {
  return(a + (b - a) * inv_logit(value));
}

vector upper_limits(vector[] a_h) {
  int nt = num_elements(a_h[1]);
  int Q = size(a_h);
  vector[nt] a_h_sums;
  vector[nt] out;
  for(k in 1:nt) {
    a_h_sums[k] = sum(a_h[1:Q, k]);
    out[k] = 1 - a_h_sums[k];
    if(out[k] <= 0) {out[k] = .00001;}
  }

  return(out);
}
vector raw_limit_to_b_h_limit(vector b_h_limit_s, vector upperLimits) {
  int nt = num_elements(upperLimits);
  vector[nt] out;
  for(k in 1:nt) {
    out[k] = a_b_scale(0, upperLimits[k], b_h_limit_s[k]);
  }
  return(out);
}
vector[] simplex_to_bh(vector[] b_h_simplex, vector b_h_limit) {
  int nt = size(b_h_simplex);
  int P = num_elements(b_h_simplex[1]);
  vector[nt] b_h[P];
  for(k in 1:nt) {
    b_h[1:P, k] = to_array_1d(b_h_simplex[k] * b_h_limit[k]);
  }
  return(b_h);
}
