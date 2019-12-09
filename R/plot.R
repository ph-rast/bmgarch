##' @title Plot 
##' @param object stan object
##' @param type Plot past expected means ("means"), or past conditional volatilyt ("cvar"), or past conditinoal correlation ("ccor")
##' @param CrI Lower and upper bound of credible interval. Default is "c( 0.025, .975)".
##' @return Plot
##' @author Philippe Rast <rast.ph@gmail.com>
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_ribbon
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 coord_cartesian
##' @export
plot.bmgarch <- function(object, type = "means", askNewPage =  TRUE, CrI = c(.025, .975)) {
    plt <- list()
    if ( type == 'means' ) {
        estimated_mean <- array(NA, dim = c(object$TS_length, object$nt))
        CIs <- array(NA, dim = c(object$TS_length, object$nt))
        for ( i in 1:object$nt ) {
            estimated_mean[,i] <- colMeans(rstan::extract(object$model_fit)[['mu']][,,i])
            ## Note: CIs will be overwritten - only stores CI's for current i
            CIs <- apply(rstan::extract(object$model_fit)[['mu']][,,i], 2, bmgarch:::.qtile, CrI  )
            df <- data.frame(mu = estimated_mean[,i], CIu = CIs[1,], CIl = CIs[2,])
            df$period <- seq_len( nrow(df) )
            plt[[i]] <- ggplot(data = df, aes(x = period, y = mu))  + labs(y = 'Conditional Means',
                                                                          x = 'Time Period',
                                                                          title = paste0(object$param, "-MGARCH"),
                                                                          subtitle = object$TS_names[[i]])
            plt[[i]] <- plt[[i]] + geom_ribbon(aes(ymin = CIl, ymax = CIu), alpha = .3 )+ geom_line()
            plot(plt[[i]])
            if ( i == 1 ) {
                devAskNewPage( ask = askNewPage )
            }
        }
    } else {
        if ( type == 'ccor' ) {
             if (object$param == 'CCC') {
                 warning('Correlation is constant for CCC - no plot generated') 
             } else {
                 ## only makes sense if model is not CCC
                 ##
                 ## Create matrix for indexing plots:
                 plot_index <- diag( object$nt )
                 ncorr <-  ( object$nt^2 - object$nt )/2
                 plot_index[upper.tri(plot_index )] <- seq_len(ncorr )
                 
                for( i in 1:(object$nt-1) ) {
                    for( j in (i+1):object$nt ) {
                        cond_corr <- apply(rstan::extract(object$model_fit)[['corH']][ , , i, j], 2, mean )
                        ci_corr <- apply(rstan::extract(object$model_fit)[['corH']][ , , i, j], 2, bmgarch:::.qtile, CrI )
                        df <- data.frame(cond_corr, lower = ci_corr[1,], upper = ci_corr[2,] )
                        df$period <- seq_len( nrow(df) )
                        df[1,] <- NA
                        plt[[ plot_index[i, j] ]] <- ggplot(data = df, aes(x = period, y = cond_corr) ) +
                            labs(y = 'Conditional Correlation',
                                 x = 'Time Period',
                                 title = paste0(object$param, "-MGARCH"),
                                 subtitle = paste( object$TS_names[[i]], object$TS_names[[j]], sep = '-') ) +
                            geom_ribbon(aes(ymin = lower, ymax = upper),
                                        alpha = .3, fill = 'red') + geom_line() +
                            coord_cartesian( ylim = c(-1, 1) )
                        
                        plot(plt[[ plot_index[i, j] ]])
                        if ( j == i+1 ) {
                            devAskNewPage( ask = askNewPage )
                        }                           
                    }
                }
             }
        } else {
            if ( type == 'cvar' ) {
                for ( i in 1:object$nt ) {
                    ## average conditional variance across iterations
                    mean_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, mean )
                    ci_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, bmgarch:::.qtile, CrI  )
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
                        devAskNewPage( ask = askNewPage )
                    }
                 }
            } 
         }
    }

    ## Return plt for use in other functions
    return_plot = list( past_data = df,
                       retro_plot = plt)
    class(return_plot) <- "plot.bmgarch"
    invisible(return_plot)
}
