if( compute_log_lik ) {
  if (distribution == 0 ) {
    for(i in 1:ahead) {
      log_lik[i] = multi_normal_lpdf( future_rts[i] | mu_forecasted[i], H_forecasted[i]);
    }
  } else if (distribution == 1 ) {
    for(i in 1:ahead) {
      log_lik[i] = multi_student_t_lpdf( future_rts[i] | nu, mu_forecasted[i], H_forecasted[i]);
    }
  }
 }
