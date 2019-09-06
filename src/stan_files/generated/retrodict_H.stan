// retrodict given distribution type
if ( distribution == 0 ){ 
  for (t in 1:T) {
    rts_out[,t] = multi_normal_rng(mu[t,], H[t,]);
    corH[t,] = cov2cor(H[t,]);
    log_lik[t] = multi_normal_lpdf(rts[t,] | mu[t,], H[t,]);
  }
 } else if ( distribution == 1 ) {
  for (t in 1:T) {
    rts_out[,t] = multi_student_t_rng(nu, mu[t,], H[t,]);
    corH[t,] = cov2cor(H[t,]);
    log_lik[t] = multi_student_t_lpdf(rts[t,] | nu, mu[t,], H[t,]);
  }
 }
