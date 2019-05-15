##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param x stan objec
##' @keywords internal
qtile = function(x){
  cis = quantile(x, c(.025, .975) )
  return(cis)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot 
##' @param object stan object
##' @param type Plot past expected means (`means`), or past conditional volatilyt (`cvar'), or past conditinoal covariance ('cov')
##' @return Plot
##' @author philippe
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 labs
plot.bmgarch = function(object, type = "means"){
    plt = list()
    if ( type == 'means' ) {
        size <- dim(rstan::extract(object$model_fit)[['rts_p']])
        ## obtain predicted values from all chains 
        forecasted_ts = array(NA, dim = c(size[1], object$nt))
        estimated_mean = array(NA, dim = c(object$TS_length, object$nt))
        CIs = array(NA, dim = c(object$TS_length, object$nt))
        for ( i in 1:object$nt ){
            forecasted_ts[,i] <- array(rstan::extract(object$model_fit)[['rts_p']][,1:object$ahead, i],
                                       dim = c(size[1], size[2]))
             estimated_mean[,i] <- colMeans(rstan::extract(object$model_fit)[['mu']][,,i])
            ## Note: CIs will be overwritten - only stores CI's for current i
            CIs = apply(rstan::extract(object$model_fit)[['mu']][,,i], 2, qtile )
            df = data.frame(mu = estimated_mean[,i], CIu = CIs[1,], CIl = CIs[2,])
            df$period = 1:nrow(df)
            plt[[i]] = ggplot(data = df, aes(x = period, y = mu))  + labs(y = 'Conditional Means',
                                                                          x = 'Time Period',
                                                                          title = paste0(object$param, "-MGARCH"),
                                                                          subtitle = object$TS_names[[i]])
            plt[[i]] = plt[[i]] + geom_ribbon(aes(ymin = CIl, ymax = CIu), alpha = .3 )+ geom_line()
            plot(plt[[i]])
            if ( i == 1){
                devAskNewPage( ask = TRUE )
            }
        }
    } else {
        if ( type == 'ccov' ) {
            ## TODO
#            dim(rstan::extract(object$model_fit)[['H']])
#            ccov_posterior =
#                rstan::extract(object$model_fit)[['H']][1, 1, 1, 1] ## [, , variance, covariance(s)...]
#            dim(ccov_posterior)
            return('To do - when nt > 2')
        } else {
            if ( type == 'cvar' ) {
                 for ( i in 1:object$nt ) {
                    mean_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, mean )
                    ci_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, qtile )
                    df = data.frame(mean_var, lower = ci_var[1,], upper = ci_var[2,] )
                    df$period = 1:nrow(df)
                    df[1,] = NA
                    plt[[i]] =  ggplot(data = df, aes(x = period, y = mean_var) ) + labs(y = 'Conditional Variance',
                                                                                         x = 'Time Period',
                                                                                         title = paste0(object$param, "-MGARCH"),
                                                                                         subtitle = object$TS_names[[i]])
                    plt[[i]] = plt[[i]] + geom_ribbon(aes(ymin = lower, ymax = upper),
                                                      alpha = .3, fill = 'red') + geom_line()
                    plot(plt[[i]])
                    if ( i == 1){
                        devAskNewPage( ask = TRUE )
                    }
                 }
             }
         }
    }
}


    
