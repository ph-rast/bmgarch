##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot 
##' @param object 
##' @return Plot
##' @author philippe
plot.bmgarch = function(object){
    size <- dim(rstan::extract(object$model_fit)[['rts_p']])
    ## obtain predicted values from all chains 
    predval = NA
    for ( i in 1:object$nt ){
        predval[i] <- array(rstan::extract(object$model_fit)[['rts_p']][,1:object$ahead, i], dim = c(size[1], size[2]))
    }
    return('yay')
}
