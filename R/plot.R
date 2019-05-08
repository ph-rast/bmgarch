##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param x 
##' @keywords internal
qtile = function(x){
  cis = quantile(x, c(.025, .975) )
  return(cis)
}


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
    forecasted_ts = array(NA, dim = c(size[1], object$nt))
    estimated_mean = array(NA, dim = c(object$TS_length, object$nt))
    for ( i in 1:object$nt ){
        forecasted_ts[,i] <- array(rstan::extract(object$model_fit)[['rts_p']][,1:object$ahead, i], dim = c(size[1], size[2]))
        estimated_mean[,i] <- colMeans(rstan::extract(object$model_fit)[['mu']][,,i])
    }
    
    CIs = apply(rstan::extract(object$model_fit)[['mu']][,,1], 2, qtile )
    df = data.frame(mu = estimated_mean, CIu = CIs[1,], CIl = CIs[2,])
    df$period = 1:nrow(df)

    plt = ggplot2::ggplot(data = df, aes(x = period, y = mu.1)) + geom_line()
    plt = plt + geom_ribbon(aes(ymin = CIl, ymax = CIu), alpha = .3, )
    return(plt)
}
