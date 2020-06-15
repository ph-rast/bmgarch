if ( distribution == 0 ) {
  rts_p[t,] = multi_normal_rng( mu_p[t,], H_p[t,]);
  if( compute_log_lik ) {
    for( i in 1:ahead ){
      log_lik[i] = multi_normal_lpdf( future_rts[i] |  mu_p[t,], H_p[t,] );
    }
  }
 } else if ( distribution == 1 ) {
  rts_p[t,] = multi_student_t_rng( nu, mu_p[t,], H_p[t,]);
  if( compute_log_lik ) {
    for( i in 1:ahead ){
      log_lik[i] = multi_student_t_lpdf( future_rts[i] | nu, mu_p[t,], H_p[t,] );
    }
  }
 }
