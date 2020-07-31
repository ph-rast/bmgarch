if ( distribution == 0 ) {
  rts_p[t,] = multi_normal_rng( mu_p[t,], H_p[t,]);
 } else if ( distribution == 1 ) {
  rts_p[t,] = multi_student_t_rng( nu, mu_p[t,], H_p[t,]);
 }
