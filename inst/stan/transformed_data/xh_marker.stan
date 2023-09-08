// Check whether xC contains a predictor or not.
array[T] matrix[nt, nt] xC_m;
int<lower = 0> xC_marker = 0;
real<lower = 0> cp;

for( t in 1:T ){
  xC_m[t] = diag_matrix( xC[t] );
  // add a variable that notes if xC is null or actually a predictor
  cp = sum( xC_m[t]' * xC_m[t] );
  if( cp != 0 )
    xC_marker = xC_marker + 1;
 }
