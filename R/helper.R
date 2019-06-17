## obtain MA piece for Forecasting
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multiply matrices in array with a vector 
##' @param MA 
##' @param theta 
##' @param mu 
##' @param rts 
##' @param i 
##' @return matrix
##' @author philippe
##' @keywords internal
.f_MA = function(MA, theta, mu, rts, i){
    res = ( -sweep( mu, 2, rts ) )
    MA[i,] = t( theta[i, ,] %*% res[i,] )
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multiply matrices in array with a vector -- generic
##' @param mat_out 
##' @param array_obj 
##' @param mat_obj 
##' @param i 
##' @return matrix
##' @author philippe
##' @keywords internal
.f_array_x_mat = function(mat_out, array_obj, mat_obj, i){
    mat_out[i,] = t( array_obj[i, ,] %*% mat_obj[i,] )
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Internal function to be used in sweep()
##' @param x Value to be squared
##' @return Squared value
##' @author philippe
##' @keywords internal
.square = function(x){
    x^2
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param x stan objec
##' @keywords internal
.cp = function(x){
    cls = length( x )
    x_a = array( x, dim = c(cls, 1) )
    x_a %*% t( x_a )
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param x stan objec
##' @keywords internal
.qtile = function(x){
  cis = quantile(x, c(.025, .975) )
  return(cis)
}
