## obtain MA piece for Forecasting
##' @title Multiply matrices in array with a vector 
##' @param MA 
##' @param theta 
##' @param mu 
##' @param rts 
##' @param i 
##' @return matrix
##' @author Philippe Rast
##' @keywords internal
.f_MA = function(MA, theta, mu, rts, i){
    res = ( -sweep( mu, 2, rts ) )
    MA[i,] = t( theta[i, ,] %*% res[i,] )
}

##' @title Multiply matrices in array with a vector -- generic
##' @param mat_out 
##' @param array_obj 
##' @param mat_obj 
##' @param i 
##' @return matrix
##' @author Philippe Rast
##' @keywords internal
.f_array_x_mat = function(mat_out, array_obj, mat_obj, i){
    mat_out[i,] = t( array_obj[i, ,] %*% mat_obj[i,] )
}

##' @title Internal function to be used in sweep()
##' @param x Value to be squared
##' @return Squared value
##' @author Philippe Rast
##' @keywords internal
.square = function(x){
    x^2
}

##' @title Internal function 
##' @param x stan objec
##' @keywords internal
.cp = function(x){
    cls = length( x )
    x_a = array( x, dim = c(cls, 1) )
    x_a %*% t( x_a )
}

##' @title Internal function to be used
##' @param x 
##' @keywords internal
##' @importFrom stats quantile
.qtile <- function(x, CrI = c(.025, .975) ) {
  cis <- quantile(x, CrI )
  return(cis)
}
