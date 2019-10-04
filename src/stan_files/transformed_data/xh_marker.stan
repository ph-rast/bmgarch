// Check whether xH contains a predictor or not.
matrix[nt, nt] xH_m[T];
int<lower = 0> xH_marker = 0;
real<lower = 0> cp;
for( t in 1:T ){
  xH_m[t] = diag_matrix( xH[t] );
  // add a variable that notes if xH is null or actually a predictor
  cp = sum( xH[t] );
  if( cp != 0)
    xH_marker = xH_marker + 1;
 }
