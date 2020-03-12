##' @title Plot method for bmgarch objects.
##' @param x bmgarch object.
##' @param type String (Default: "mean"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
##' @param askNewPage askNewPage Logical (Default: True). Whether to ask for new plotting page.
##' @param CrI CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @return List of ggplot objects (one per time series).
##' @author Stephen R. Martin
##' @import ggplot2
##' @importFrom graphics plot
##' @importFrom grDevices devAskNewPage
##' @export
plot.bmgarch <- function(x, type = "mean", askNewPage = TRUE, CrI = c(.025, .975)) {
    x.fitted <- fitted(x, CrI = CrI, digits = 4)
    x.observed <- x$RTS_full
    nt <- x$nt
    TS_names <- lapply(x.fitted$backcast, function(x) {dimnames(x)[[3]]})
    TS_length <- x$TS_length

    df <- as.data.frame(x.fitted)

    # rename % cols to L and U.
    LU <- paste0(CrI * 100, "%")
    colnames(df)[colnames(df) %in% LU] <- c("L", "U")

    plt <- list()

    # Subset by correct param type.
    if(!(type %in% c("mean","var","cor"))) {
        stop("'type' must be 'mean', 'var', or 'cor'.")
    }
    if(type == "cor" & x$param == "CCC") {
        stop("CCC does not model correlations over time.")
    }

    df <- df[df$param == type, ]

    for(i in TS_names[[type]]) {
        df.i <- df[df$TS == i,]

        plt[[i]] <- ggplot(data = df.i, aes(x = period, y = mean, ymin = L, ymax = U)) +
            geom_line() +
            geom_ribbon(alpha = .3)
        l <- labs(x = "Time Period",
                  y = switch(type,
                            mean = "Conditional Means",
                            var = "Conditional Variances",
                            cor = "Conditional Correlations",
                            NULL),
                  title = paste0(x$param, "-MGARCH"),
                  subtitle = i)
        plt[[i]] <- plt[[i]] + l
        print(plt[[i]])
        if(which(i %in% TS_names[[type]]) == 1) {
            devAskNewPage(ask = askNewPage)
        }
    }

    return(invisible(plt))
}
##' @title Plot method for forecast.bmgarch objects.
##' @param x forecast.bmgarch object. See \code{\link{forecast.bmgarch}}.
##' @param type String (Default: "mean"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
##' @param askNewPage Logical (Default: True). Whether to ask for new plotting page.
##' @param last_t Integer (Default: 100). For plotting only. Only show \code{last_t} observations in plot.
##' @return List of ggplot objects (one per time series).
##' @author Stephen R. Martin
##' @import ggplot2
##' @importFrom graphics plot
##' @importFrom grDevices devAskNewPage
##' @export
plot.forecast.bmgarch <- function(x, type = "mean", askNewPage = TRUE, last_t = 100) {
    nt <- x$meta$nt
    TS_names <- lapply(x$forecast, function(x) {dimnames(x)[[3]]})
    x.observed <- as.data.frame(x$meta$RTS_full)
    x.observed$period <- seq_len(x$meta$TS_length)

    df <- as.data.frame(x, backcast = TRUE)

    # rename % cols to L and U.
    LU <- paste0(x$meta$CrI * 100, "%")
    colnames(df)[colnames(df) %in% LU] <- c("L", "U")

    # Subset by correct param type.
    if(!(type %in% c("mean","var","cor"))) {
        stop("'type' must be 'mean', 'var', or 'cor'.")
    }
    if(type == "cor" & x$meta$param == "CCC") {
        stop("CCC does not model correlations over time.")
    }
    df <- df[df$param == type, ]
    df$type <- ifelse(df$type == "backcast", "Backcast", "Forecast")

    plt <- list()
    TS_length <- max(df$period)

    for(i in TS_names[[type]]) {
        df.i <- df[df$TS == i,]

        plt[[i]] <- ggplot(data = df.i, aes(x = period, y = mean, ymin = L, ymax = U, color = type, fill = type)) +
            geom_line() +
            geom_ribbon(alpha = .3)
        l <- labs(x = "Time Period",
                  y = switch(type,
                            mean = "Conditional Means",
                            var = "Conditional Variances",
                            cor = "Conditional Correlations",
                            NULL),
                  title = paste0(x$meta$param, "-MGARCH"),
                  subtitle = i,
                  color = "Type",
                  fill = "Type")
        plt[[i]] <- plt[[i]] + l
        plt[[i]] <- plt[[i]] + coord_cartesian(xlim = c(TS_length - last_t, TS_length))

        print(plt[[i]])
        if(which(i %in% TS_names[[type]]) == 1) {
            devAskNewPage(ask = askNewPage)
        }
    }

    return(invisible(plt))
}


## @title Plot method for bmgarch objects.
## @param x bmgarch object.
## @param type String (Default: "mean"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
## @param askNewPage Logical (Default: True). Whether to ask for new plotting page.
## @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
## @return plot.bmgarch (Invisibly). List of plots.
## @author Philippe Rast
## @importFrom ggplot2 ggplot
## @importFrom ggplot2 geom_line
## @importFrom ggplot2 geom_ribbon
## @importFrom ggplot2 aes
## @importFrom ggplot2 labs
## @importFrom ggplot2 coord_cartesian
## @importFrom graphics plot
## @importFrom grDevices devAskNewPage
## @export
## plot.bmgarch.old <- function(x, type = "mean", askNewPage =  TRUE, CrI = c(.025, .975)) {
##     object <- x
##     plt <- list()
##     if ( type == 'mean' ) {
##         estimated_mean <- array(NA, dim = c(object$TS_length, object$nt))
##         CIs <- array(NA, dim = c(object$TS_length, object$nt))
##         for ( i in 1:object$nt ) {
##             estimated_mean[,i] <- colMeans(rstan::extract(object$model_fit)[['mu']][,,i])
##             ## Note: CIs will be overwritten - only stores CI's for current i
##             CIs <- apply(rstan::extract(object$model_fit)[['mu']][,,i], 2, .qtile, CrI  )
##             df <- data.frame(mu = estimated_mean[,i], CIu = CIs[1,], CIl = CIs[2,])
##             df$period <- seq_len( nrow(df) )
##             plt[[i]] <- ggplot(data = df, aes(x = period, y = mu))  + labs(y = 'Conditional Means',
##                                                                           x = 'Time Period',
##                                                                           title = paste0(object$param, "-MGARCH"),
##                                                                           subtitle = object$TS_names[[i]])
##             plt[[i]] <- plt[[i]] + geom_ribbon(aes(ymin = CIl, ymax = CIu), alpha = .3 )+ geom_line()
##             plot(plt[[i]])
##             if ( i == 1 ) {
##                 devAskNewPage( ask = askNewPage )
##             }
##         }
##     } else {
##         if ( type == 'cor' ) {
##              if (object$param == 'CCC') {
##                  warning('Correlation is constant for CCC - no plot generated') 
##              } else {
##                  ## only makes sense if model is not CCC
##                  ##
##                  ## Create matrix for indexing plots:
##                  plot_index <- diag( object$nt )
##                  ncorr <-  ( object$nt^2 - object$nt )/2
##                  plot_index[upper.tri(plot_index )] <- seq_len(ncorr )
                 
##                 for( i in 1:(object$nt-1) ) {
##                     for( j in (i+1):object$nt ) {
##                         cond_corr <- apply(rstan::extract(object$model_fit)[['corH']][ , , i, j], 2, mean )
##                         ci_corr <- apply(rstan::extract(object$model_fit)[['corH']][ , , i, j], 2, .qtile, CrI )
##                         df <- data.frame(cond_corr, lower = ci_corr[1,], upper = ci_corr[2,] )
##                         df$period <- seq_len( nrow(df) )
##                         df[1,] <- NA
##                         plt[[ plot_index[i, j] ]] <- ggplot(data = df, aes(x = period, y = cond_corr) ) +
##                             labs(y = 'Conditional Correlation',
##                                  x = 'Time Period',
##                                  title = paste0(object$param, "-MGARCH"),
##                                  subtitle = paste( object$TS_names[[i]], object$TS_names[[j]], sep = '-') ) +
##                             geom_ribbon(aes(ymin = lower, ymax = upper),
##                                         alpha = .3, fill = 'red') + geom_line() +
##                             coord_cartesian( ylim = c(-1, 1) )
                        
##                         plot(plt[[ plot_index[i, j] ]])
##                         if ( j == i+1 ) {
##                             devAskNewPage( ask = askNewPage )
##                         }                           
##                     }
##                 }
##              }
##         } else {
##             if ( type == 'var' ) {
##                 for ( i in 1:object$nt ) {
##                     ## average conditional variance across iterations
##                     mean_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, mean )
##                     ci_var = apply(rstan::extract(object$model_fit)[['H']][ , , i, i], 2, .qtile, CrI  )
##                     df = data.frame(mean_var, lower = ci_var[1,], upper = ci_var[2,] )
##                     df$period = 1:nrow(df)
##                     df[1,] = NA
##                     plt[[i]] =  ggplot(data = df, aes(x = period, y = mean_var) ) + labs(y = 'Conditional Variance',
##                                                                                          x = 'Time Period',
##                                                                                          title = paste0(object$param, "-MGARCH"),
##                                                                                          subtitle = object$TS_names[[i]])
##                     plt[[i]] = plt[[i]] + geom_ribbon(aes(ymin = lower, ymax = upper),
##                                                       alpha = .3, fill = 'red') + geom_line()
##                     plot(plt[[i]])
##                     if ( i == 1){
##                         devAskNewPage( ask = askNewPage )
##                     }
##                  }
##             } 
##          }
##     }

##     ## Return plt for use in other functions
##     return_plot = list( past_data = df,
##                        retro_plot = plt)
##     class(return_plot) <- "plot.bmgarch"
##     invisible(return_plot)
## }
