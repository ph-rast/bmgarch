if( meanstructure == 0 ){
  // Constant
  mu[t, ] = phi0;
 } else if( meanstructure == 1 ){
  // arma11
  mu[t, ] = phi0 + phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
 } else if( meanstructure == 2 ){
  // VAR1
  mu[t, ] = phi0 + phi * rts[t-1, ];
 }
