if( meanstructure == 0 ){
  mu[t, ] = phi0;
 } else if( meanstructure == 1 ){
  mu[t, ] = phi0 + phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
 }

