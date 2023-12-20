/*
  Computes jacobian for the a_b_scale transform.
  @param a Lower limit.
  @param b Upper limit.
  @param value The unconstrained value to be transformed.
  @return log jacobian contribution.
 */
real a_b_scale_jacobian(real a, real b, real value) {
  real invlogit_value = inv_logit(value);
  real out = log(b - a) + log(invlogit_value) + log1m(invlogit_value);
  return(out);
}

/*
  Takes unconstrained R value, transforms to be between a and b.
  @param a Lower limit.
  @param b Upper limit.
  @param value The unconstrained value to transform
  @return real value between a and b.
 */
real a_b_scale(real a, real b, real value) {
  return(a + (b - a) * inv_logit(value));
}

/*
  Compute the upper limit of b_h from a_h, given that sum(a_h) + sum(b_h) < 1
  This works for b_h given a_h, or a_h given b_h. They yield identical models, assuming a dual constraint, such that 0 < sum(a_h) < 1 and 0 < sum(b_h) < 1 - sum(a_h), or vice-versa.
  @param a_h Array of vectors [nt, Q].
  @return vector[nt] Upper limits such that sum(b_h{k}) < UpperLimit{k} = 1 - sum(a_h{k})
 */
vector upper_limits(array[] vector a_h) {
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

/*
  Transform unconstrained sums to (0, upperLimits).
  @param b_h_sum_s A vector [nt] of unconstrained sums.
  @param upperLimits A vector [nt] of the upper limits.
  @return vector[nt] of values between [0, upperLimits].
 */
vector raw_sum_to_b_h_sum(vector b_h_sum_s, vector upperLimits) {
  int nt = num_elements(upperLimits);
  vector[nt] out;
  for(k in 1:nt) {
    out[k] = a_b_scale(0, upperLimits[k], b_h_sum_s[k]);
  }
  return(out);
}
/*
  Transform the array of b_h simplexes to b_h.
  
  The dual constraint is such that 0 < sum(a_h) + sum(b_h) < 1, for each timeseries.
  We assume sum(a_h) is [0,1], then sum(b_h) is [0, 1 - sum(a_h)].
  The sums are estimated for a_h and b_h.
  Then the a_h and b_h simplexes are scaled by the sums to produce a_h and b_h.

  @param b_h_simplex Array [nt] of simplexes [P].
  @param b_h_sum vector [nt] of the value to which the b_h's should be summed.
  @return An array [P] of vectors [nt].
 */
array[] vector simplex_to_bh(array[] vector b_h_simplex, vector b_h_sum) {
  int nt = size(b_h_simplex);
  int P = num_elements(b_h_simplex[1]);
  array[P] vector[nt] b_h;
  for(k in 1:nt) {
    b_h[1:P, k] = to_array_1d(b_h_simplex[k] * b_h_sum[k]);
  }
  return(b_h);
}
