##' Estimates forecasted means, variances, and correlations from a fitted bmgarch model.
##' @title Forecast method for bmgarch objects.
##' @param object bmgarch object.
##' @param ahead Integer (Default: 1). Periods to be forecasted ahead.
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c. Used in a log-linear model on the constant variance terms. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param seed Integer (Optional). Specify seed for \code{\link[rstan]{sampling}}.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @return forecast.bmgarch object. List containing \code{forecast}, \code{backcast}, and \code{meta}data. See \code{\link{fitted.bmgarch}} for information on \code{backcast}.
##' \code{forecast} is a list of four components. Three components (mean, var, cor) contain forecasted values in arrays. mean and var are structured as \code{[N, 7, TS]} arrays, where \code{N} is length of the time series, and \code{TS} is the number of time series. \code{cor} is structured as \code{[N, 7, TS_TS]}, where \code{TS_TS} represents a (lower triangular) correlation. All arguments are named. E.g., \code{forecasted$forecast$cor[, , "tsB_tsA"]} would the forecasted values for the correlation between the two time series, "tsB" and "tsA".
##' @author Stephen R. Martin
##' @importFrom forecast forecast
##' @export forecast
##' @export
forecast.bmgarch <- function(object, ahead = 1, xC = NULL, CrI = c(.025, .975), seed = NA, digits = 2) {

    # Define a 0 array for stan.
    if(is.null(xC)) {
        xC <- array(0, dim = c(ahead, object$nt))
    }

    # Get stan data
    standat <- list(T = object$TS_length,
                     nt = object$nt,
                     rts = cbind(object$RTS_full),
                     xC = object$xC,
                     Q =  object$mgarchQ,
                     P =  object$mgarchP,
                     ahead =  ahead, 
                     meanstructure =  object$meanstructure,
                     distribution =  object$num_dist,
                     xC_p =  xC)

    gqs_model <- switch(object$param,
                        DCC = stanmodels$forecastDCC,
                        CCC = stanmodels$forecastCCC,
                        BEKK = stanmodels$forecastBEKK,
                        pdBEKK = stanmodels$forecastBEKK,
                        NULL)

    if(is.null(gqs_model)) {
        stop("bmgarch object 'param' does not match a supported model. ",
             object$param, "is not one in ", paste0(supported_models, collapse = ", "), ".")
    }

    backcast <- max(object$mgarchP, object$mgarchQ)
    nt <- object$nt
    cast_start <- (object$TS_length - backcast + 1)
    forecast_start <- (object$TS_length + 1)
    forecast_end <- (object$TS_length + ahead)

    # TODO: Limit pars to only what is needed (H_p, R/R_p, rts_p, mu_p)
    forecasted <- rstan::gqs(gqs_model,
                             draws = as.matrix(object$model_fit),
                             data = standat,
                             seed = seed)

    # Get stan summaries for forecasted values.
    f.mean <- .get_stan_summary(forecasted, "rts_p", CrI)
    f.var <- .get_stan_summary(forecasted, "H_p", CrI)
    f.cor <- NA
    if(object$param != "CCC") {
        f.cor <- .get_stan_summary(forecasted, "R_p", CrI)
    }

    # Restructure to array
    ## f.mean
    stan_sum_cols <- colnames(f.mean)
    f.mean <- array(f.mean, dim = c(nt, backcast + ahead , ncol(f.mean)))
    f.mean <- aperm(f.mean, c(2,3,1))
    dimnames(f.mean) <- list(period = cast_start:forecast_end, stan_sum_cols, TS = object$TS_names)

    ## f.var
    ### Pull out indices for [period, a, a]
    f.var.indices <- grep("H_p\\[[[:digit:]]+,([[:digit:]]+),\\1]", rownames(f.var), value = TRUE)
    f.var <- f.var[f.var.indices,]
    f.var <- array(f.var, dim = c(nt, backcast + ahead, ncol(f.var)))
    f.var <- aperm(f.var, c(2, 3, 1))
    dimnames(f.var) <- list(period = cast_start:forecast_end, stan_sum_cols, TS = object$TS_names)

    ## f.cor
    if(object$param != "CCC") {
        # Lower-triangular indices
        f.cor.indices.L <- which(lower.tri(matrix(0, nt, nt)), arr.ind = TRUE)
        # Labels mapping to TS names
        f.cor.indices.L.labels <- paste0(object$TS_names[f.cor.indices.L[,1]], "_", object$TS_names[f.cor.indices.L[,2]])
        # Indices as "a,b"
        f.cor.indices.L.char <- paste0(f.cor.indices.L[,1], ",", f.cor.indices.L[,2])
        # Indicices as "[period,a,b]"
        f.cor.indices.L.all <- paste0("R_p[",1:(backcast + ahead), ",", rep(f.cor.indices.L.char, each = (backcast + ahead)),"]")
        # Get only these elements.
        f.cor <- f.cor[f.cor.indices.L.all,]
        f.cor <- array(f.cor, dim = c(backcast + ahead, length(f.cor.indices.L.char), ncol(f.cor)))
        f.cor <- aperm(f.cor, c(1, 3, 2))
        dimnames(f.cor) <- list(period = cast_start:forecast_end, stan_sum_cols, TS = f.cor.indices.L.labels)
    }

    # Remove backcasts from forecasts.
    f.mean <- f.mean[-c(1:backcast), , , drop = FALSE]
    f.var <- f.var[-c(1:backcast), , , drop = FALSE]
    if(object$param != "CCC") {
        f.cor <- f.cor[-c(1:backcast), , , drop = FALSE]
    }

    out <- list()
    out$forecast$mean <- f.mean
    out$forecast$var <- f.var
    out$forecast$cor <- f.cor
    out$forecast$meta <- list(xC = xC, TS_length = ahead)

    metaNames <- c("param", "distribution", "num_dist", "nt", "TS_length", "TS_names", "RTS_full", "mgarchQ", "mgarchP", "xC", "meanstructure")
    meta <- with(object, mget(metaNames))
    out$meta <- meta
    out$meta$digits <- digits

    out$backcast <- fitted.bmgarch(object, CrI)$backcast

    class(out) <- "forecast.bmgarch"
    return(out)

}
##' Extracts the model-predicted means, variances, and correlations for the fitted data.
##'
##' Whereas \code{\link{forecast.bmgarch}} computes the \emph{forecasted} values for future time periods, \code{fitted.bmgarch} computes the \emph{backcasted} (model-predicted) values for the observed time periods.
##' @title Fitted (backcasting) method for bmgarch objects.
##' @param object bmgarch object.
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param ... Not used.
##' @return 
##' @author Stephen R. Martin
##' @importFrom stats fitted
##' @export
fitted.bmgarch <- function(object, CrI = c(.025, .975), ...) {
    nt <- object$nt
    TS_length <- object$TS_length

    b.mean <- .get_stan_summary(object$model_fit, "mu", CrI)
    b.var <- .get_stan_summary(object$model_fit, "H", CrI)
    b.cor <- NA
    if(object$param != "CCC") {
        ## In forecast stan, R is the same as corH.
        ## In DCC, corH is the same as R also.
        ## Not worth changing the stan files to be consistent, since structures are the same.
        b.cor <- .get_stan_summary(object$model_fit, "corH", CrI)
    }

    # Restructure
    ## b.mean
    stan_sum_cols <- colnames(b.mean)
    b.mean <- array(b.mean, dim = c(nt, TS_length, ncol(b.mean)))
    b.mean <- aperm(b.mean, c(2,3,1))
    dimnames(b.mean) <- list(period = 1:TS_length, stan_sum_cols, TS = object$TS_names)

    ## b.var
    b.var.indices <- grep("H\\[[[:digit:]]+,([[:digit:]]+),\\1]", rownames(b.var), value = TRUE)
    b.var <- b.var[b.var.indices,]
    b.var <- array(b.var, dim = c(nt, TS_length, ncol(b.var)))
    b.var <- aperm(b.var, c(2, 3, 1))
    dimnames(b.var) <- list(period = 1:TS_length, stan_sum_cols, TS = object$TS_names)

    ## b.cor
    if(object$param != "CCC") {
        # Lower-triangular indices
        b.cor.indices.L <- which(lower.tri(matrix(0, nt, nt)), arr.ind = TRUE)
        # Labels mapping to TS names
        b.cor.indices.L.labels <- paste0(object$TS_names[b.cor.indices.L[,1]], "_", object$TS_names[b.cor.indices.L[,2]])
        # Indices as "a,b"
        b.cor.indices.L.char <- paste0(b.cor.indices.L[,1], ",", b.cor.indices.L[,2])
        # Indicices as "[period,a,b]"
        b.cor.indices.L.all <- paste0("corH[",1:object$TS_length, ",", rep(b.cor.indices.L.char, each = object$TS_length),"]")
        # Get only these elements.
        b.cor <- b.cor[b.cor.indices.L.all,]
        b.cor <- array(b.cor, dim = c(object$TS_length, length(b.cor.indices.L.char), ncol(b.cor)))
        b.cor <- aperm(b.cor, c(1, 3, 2))
        dimnames(b.cor) <- list(period = 1:object$TS_length, stan_sum_cols, TS = b.cor.indices.L.labels)
    }

    out <- list()
    out$backcast$mean <- b.mean
    out$backcast$var <- b.var
    out$backcast$cor <- b.cor

    metaNames <- c("param", "distribution", "num_dist", "nt", "TS_length", "TS_names", "RTS_full", "mgarchQ", "mgarchP", "xC", "meanstructure")
    meta <- with(object, mget(metaNames))
    out$meta <- meta

    class(out) <- "fitted.bmgarch"
    return(out)
}

## Forecasts the (conditional) means, conditional variances, and conditional correlations.
##
## \code{forecast} takes a fitted \code{bmgarch} object, and predicts the next \code{ahead} means and variances, with uncertainty.
## Time-varying predictors can be included in \code{xC}.
## @title Forecast method for bmgarch objects.
## @param object bmgarch object. The fitted model used for forecasting.
## @param ahead Integer (Default: 1). Periods to be forecasted ahead.
## @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c. Used in a log-linear model on the constant variance terms. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
## @param type String (Default: "var"). Whether to plot conditional means ("mean"), variance ("var"), or correlations ("cor"). 
## @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
## Possible values are .025, .05, .10, .50, .90, .95, and .975
## @param plot Logical (Default: TRUE). Should forecasts be plotted.
## @param last_t Integer (Default: 100). For plotting only. Only show \code{last_t} observations in plot.
## @param seed Integer (Optional). Specify seed for \code{\link[rstan]{sampling}}.
## @return Forecasted conditional means, variances and correlations. Forecasts are shown in plots (unless plot = FALSE).
## Conditional variance and correlation forecasts are appended to the corresponding plots with the estimates over the given period.
## Mean forecasts are appended to the observed values if no meanstructure is estimated (default), otherwise they will be appended
## to the esimated means. 
## @author Philippe Rast, Stephen R. Martin
## @importFrom ggplot2 geom_point
## @importFrom ggplot2 geom_errorbar
## @importFrom ggplot2 aes_string
## @importFrom forecast forecast
## @export forecast
## @aliases forecast
## @export
#forecast.bmgarch.old <- function(object,
#                     ahead =  1,
#                     type =  "var",
#                     xC =  NULL,
#                     CrI =  c(.025, .975),
#                     seed = NA,
#                     plot =  TRUE,
#                     last_t =  100) {
#
#    ## Define a 0 array as the stan models need some non-null matrix
#    if( is.null( xC ) ) xC <- array( rep(0, object$nt),  dim = c(ahead,  object$nt ) )
#
#    ## obtain data to be passed to gqs
#    standat <-  list(T = object$TS_length,
#                     nt = object$nt,
#                     rts = cbind(object$RTS_full),
#                     xC = object$xC,
#                     Q =  object$mgarchQ,
#                     P =  object$mgarchP,
#                     ahead =  ahead, 
#                     meanstructure =  object$meanstructure,
#                     distribution =  object$num_dist,
#                     xC_p =  xC)
#
#    gqs_model <- switch(object$param,
#                        DCC = stanmodels$forecastDCC,
#                        CCC = stanmodels$forecastCCC,
#                        BEKK = stanmodels$forecastBEKK,
#                        pdBEKK = stanmodels$forecastBEKK,
#                        NULL)
#    if(is.null(gqs_model)) {
#        stop("bmgarch object 'param' does not match a supported model. ",
#             object$param, "is not one in ", paste0(supported_models, collapse = ", "), ".")
#    }
#
#    ## Call forecasting functions
#    forecasted <- rstan::gqs(gqs_model,
#                             draws = as.matrix(object$model_fit),
#                             data = standat,
#                             seed = seed)
#    
#    ## We only need rts_p and H_p and R_p
#    ## note that rts_p and H_p include past observations (given by P or Q)
#    ## and forecasted elements given (given by ahead)
#
#    ## Global variables
#    ## Separate estimate (backcast) from forecast, given dimension of P, Q, and ahead
#    backcast <-  max(object$mgarchP, object$mgarchQ)
#    nt <- object$nt
#
#    ## ###########################
#    ## Plot and print forecasts ##
#    ## ###########################
#
#    ## init plt
#    plt <- list()
#    
#    if( type == "mean" ) {
#        ## Two versions: When no meanstructure present (default), predicted means
#        ## will be plotted with the observed values.
#        ## If meanstructure is present ( != 0 ), forecast is added to plot with
#        ## estimated means and CrI's
#
#        ## collect rts elements
#        rts_pred_mn <-  rstan::summary(forecasted, pars =  "rts_p",
#                                       probs =  .5 )$summary[, 4]        
#        
#        rts_pred_Lcri <-  rstan::summary(forecasted, pars =  "rts_p",
#                                         probs =  CrI[1] )$summary[, 4]
#        rts_pred_Ucri <-  rstan::summary(forecasted, pars =  "rts_p",
#                                         probs =  CrI[2] )$summary[, 4]
#        
#        ## Dummy code periods: observed vs forecast
#        rts_period <- rep( seq_len( length( rts_pred_mn )/nt ), each = nt )
#        rts_frcst <- ifelse( rts_period < max(rts_period) - (ahead -1 ), 0, 1)
#        
#        rts_mn <- cbind(rts_pred_mn, rts_period, rts_frcst)
#        rts_Lcri <- cbind(rts_pred_Lcri, rts_period, rts_frcst)
#        rts_Ucri <- cbind(rts_pred_Ucri, rts_period, rts_frcst)
#
#        ## ##########
#        ## Means:  ##
#        ## ##########
#
#        ## select only forecasted
#        rts_mn_forecasted <- rts_mn[rts_mn[, 3] != 0, ]
#        rts_Lcri_forecasted <- rts_Lcri[ rts_Lcri[, 3] != 0, ]
#        rts_Ucri_forecasted <- rts_Ucri[ rts_Ucri[, 3] != 0, ]
#
#        ## Write into new object
#        forecasted_data_rts <- cbind(rts_mn_forecasted[, 1],
#                                     rts_Lcri_forecasted[, 1],
#                                     rts_Ucri_forecasted[, 1])
#
#        colnames(forecasted_data_rts) <- c('rts_forecasted',
#                                           paste0('CrI_', CrI[1]),
#                                           paste0('CrI_', CrI[2]))
#
#        forecasted_rts <- list()
#
#        for(i in 1:nt ) {
#            forecasted_rts[[i]] <- forecasted_data_rts[seq(1, ahead * nt , by = nt)+(i-1), ]            
#            names(forecasted_rts[[i]])[1] <- paste(colnames(standat$rts)[i],
#                                                      names(forecasted_rts[[i]])[1], sep = "_")
#        }
#
#        
#        ## Plots:
#        if( plot ) {
#            if ( object$meanstructure != 0 ) {                
#                ## Use conditional means plot of bmgarch::plot
#                retro <- plot(object, type = 'mean', askNewPage = FALSE, CrI = CrI )
#                
#                for( i in seq_len(nt) ) {
#                    df <- array( NA, dim = c(standat$T, 3) )
#                    df <- data.frame( rbind(df, forecasted_rts[[i]]) )
#                    df$period <- seq_len( dim(df )[1] )
#                    names(df )[1] <- object$TS_names[i]
#                    
#                    plt[[i]] <- retro$retro_plot[[i]] +
#                        geom_line(data = df, aes_string( x = "period", y = names(df)[1] ) ) +                
#                        geom_errorbar(data = df,  aes_string(x = "period", ymin = names(df)[2], ymax = names(df)[3]),
#                                      inherit.aes = FALSE , alpha =  .3)+
#                        geom_point(data = df, mapping = aes_string(x = "period", y = names(df)[1] ),
#                                   inherit.aes = FALSE)+
#                        coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
#                    
#                    plot(plt[[i]])
#                    if ( i == 1) {
#                        devAskNewPage( ask = TRUE )
#                    }
#                }
#            } else if ( object$meanstructure == 0 ) {
#                ## Plot on observed data
#                rts_data <- object$RTS_full
#                
#                for( i in seq_len(nt) ) {
#                    df <- array( NA, dim = c(standat$T, 3) )
#                    df[, 1] <- rts_data[, i]
#                    df <- data.frame( rbind(df, forecasted_rts[[i]]) )
#                    df$Type <- factor( c( rep("Observed", standat$T),
#                                         rep("Forecasted",  ahead)) )
#                    
#                    df$period <- seq_len( dim(df )[1] )
#                    names(df )[1] <- object$TS_names[i]
#                    
#                    plt[[i]] <- ggplot(data = df, aes_string(x = "period",
#                                                             y = names(df)[1],
#                                                             color =  "Type") ) +
#                        geom_line() +
#                        geom_errorbar( aes_string(ymin = names(df)[2], ymax = names(df)[3]),
#                                      alpha =  .3)+
#                        geom_point(data = df[df$Type == "Forecasted", ], mapping = aes_string(x = "period", y = names(df)[1] ))+
#                        coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
#                    
#                    plot(plt[[i]])
#                    if ( i == 1) {
#                        devAskNewPage( ask = TRUE )
#                    }
#                }  
#            }
#            
#        }
#        
#        return( forecasted_rts )
#        
#    } else if(type == "var" ) {
#        ## ############
#        ## Variances ##
#        ## ############
#        ## Print variances
#        ## obtain predicted variances:
#
#        ## collect H elements
#        H_pred_mn <-  rstan::summary(forecasted, pars =  "H_p",
#                                     probs =  c(.5) )$summary[, 4]
#        H_pred_Lcri <-  rstan::summary(forecasted, pars =  "H_p",
#                                       probs =  CrI[1] )$summary[, 4]
#        H_pred_Ucri <-  rstan::summary(forecasted, pars =  "H_p",
#                                       probs =  CrI[2] )$summary[, 4]
#        ## Dummy code periods: observed vs forecast
#        H_period <- rep( seq_len( length( H_pred_mn )/nt^2 ), each = nt^2 )
#        H_frcst <- ifelse( H_period < max(H_period) - (ahead -1 ), 0, 1)
#        
#        H_mn <- cbind(H_pred_mn, H_period, H_frcst)
#        H_Lcri <- cbind(H_pred_Lcri, H_period, H_frcst)
#        H_Ucri <- cbind(H_pred_Ucri, H_period, H_frcst)
#        
#        ## select only forecasted
#        H_mn_forecasted <- H_mn[H_mn[, 3] != 0, ]
#        H_Lcri_forecasted <- H_Lcri[ H_Lcri[, 3] != 0, ]
#        H_Ucri_forecasted <- H_Ucri[ H_Ucri[, 3] != 0, ]
#        
#        ## find variances
#        
#        splitted <- strsplit( rownames(H_mn_forecasted), split = ",|]")
#        row <- as.numeric( t( sapply( splitted, "[", 2:3) )[,1])
#        col <- as.numeric( t( sapply( splitted, "[", 2:3) )[,2])
#        
#        ## obtain position of variances:
#        var_pos <- which( row == col )
#
#        ## Write into new object
#        forecasted_data_H <- cbind(H_mn_forecasted[var_pos, 1],
#                                   H_Lcri_forecasted[var_pos, 1],
#                                   H_Ucri_forecasted[var_pos, 1])
#
#        colnames(forecasted_data_H) <- c('H_forecasted',
#                                         paste0('CrI_', CrI[1]),
#                                         paste0('CrI_', CrI[2]))
#
#        forecasted_H <- list()
#        
#        for(i in 1:nt ) {
#            forecasted_H[[i]] <- forecasted_data_H[seq(1, ahead * nt , by = nt)+(i-1), ]            
#            names(forecasted_H[[i]])[1] <- paste(colnames(standat$rts)[i],
#                                                    names(forecasted_H[[i]])[1], sep = "_")
#        }
#        
#        ## Plots
#        if( plot ) {
#            ## Use conditional variance plot of bmgarch::plot
#            retro <- plot(object, type = 'var', askNewPage = FALSE, CrI = CrI)
#
#            ## Loop through nt's
#            for(i in seq_len(nt) ) {
#                df <- array( NA, dim = c(standat$T, 3) )
#                df <- data.frame( rbind(df, forecasted_H[[i]]) )
#                df$period <- seq_len( dim(df )[1] )
#                df$Type <- factor( c( rep("Estimated", standat$T),
#                                     rep("Forecasted",  ahead)) )
#
#                
#                plt[[i]] <- retro$retro_plot[[i]] +
#                    geom_line(data = df, aes_string( x = "period", y = names(df)[1],  color =  "Type" ) ) +
#                    geom_errorbar(data = df,  aes_string(x = "period", ymin = names(df)[2], ymax = names(df)[3]),
#                                  inherit.aes = FALSE , alpha =  .3)+
#                    geom_point(data = df, mapping = aes_string(x = "period",
#                                                               y = names(df)[1],
#                                                               color =  "Type" ), inherit.aes = FALSE)+
#                    coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
#                
#                plot(plt[[i]])
#                if ( i == 1 ) {
#                    devAskNewPage( ask = TRUE)
#                }
#            }
#        }
#        ## print
#        return( forecasted_H )
#    } else if(type == "cor" ) {
#        ## ############
#        ## Correlations  ##
#        ## ############
#
#        ## collect R elements
#        R_pred_mn <-  rstan::summary(forecasted, pars =  "R_p",
#                                     probs =  c(.5) )$summary[, 4]
#        R_pred_Lcri <-  rstan::summary(forecasted, pars =  "R_p",
#                                       probs =  CrI[1] )$summary[, 4]
#        R_pred_Ucri <-  rstan::summary(forecasted, pars =  "R_p",
#                                       probs =  CrI[2] )$summary[, 4]
#
#        ## Dummy code periods: observed vs forecast
#        R_period <- rep( seq_len( length( R_pred_mn )/nt^2 ), each = nt^2 )
#        R_frcst <- ifelse( R_period < max(R_period) - (ahead -1 ), 0, 1)
#        
#        R_mn <- cbind(R_pred_mn, R_period, R_frcst)
#        R_Lcri <- cbind(R_pred_Lcri, R_period, R_frcst)
#        R_Ucri <- cbind(R_pred_Ucri, R_period, R_frcst)
#
#        ## select only forecasted
#        R_mn_forecasted <- R_mn[R_mn[, 3] != 0, ]
#        R_Lcri_forecasted <- R_Lcri[ R_Lcri[, 3] != 0, ]
#        R_Ucri_forecasted <- R_Ucri[ R_Ucri[, 3] != 0, ]
#        
#        ## find correlations
#        splitted <- strsplit( rownames(R_mn_forecasted), split = ",|]")
#        row <- as.numeric( t( sapply( splitted, "[", 2:3) )[,1])
#        col <- as.numeric( t( sapply( splitted, "[", 2:3) )[,2])
#
#        ## obtain position of correlations:
#        cor_pos <- which( row != col  & row < col )
#        
#        row_col <- cbind(row, col )[cor_pos, ]
#
#        ## Write int new object
#        forecasted_data_R <- cbind(R_mn_forecasted[cor_pos, 1],
#                                   R_Lcri_forecasted[cor_pos, 1],
#                                   R_Ucri_forecasted[cor_pos, 1])
#
#        colnames(forecasted_data_R) <- c('R_frcsted',
#                                         paste0('CrI_', CrI[1]),
#                                         paste0('CrI_', CrI[2]))
#
#        forecasted_R <- list()
#        corrs <- (nt^2 - nt)/2
#        for(i in 1:corrs  ) {
#            forecasted_R[[i]] <- forecasted_data_R[seq(1, ahead * corrs , by = corrs)+(i-1), ]
#            if(corrs == 1 ) label_pos <-  1 else {
#            label_pos <- row_col[seq(1, ahead * corrs , by = corrs)+(i-1), ][1, ]}
#            names(forecasted_R[[i]])[1] <- paste( paste(abbreviate( colnames(standat$rts)[label_pos[1]]),
#                                                        abbreviate( colnames(standat$rts)[label_pos[2]]),
#                                                           sep =  "_"),
#                                                    names(forecasted_R[[i]])[1], sep = "_")
#        }
#        
#        ## Plots
#        if( plot ) {
#            ## Use conditional variance plot of bmgarch::plot
#            retro <- plot(object, type = 'cor', askNewPage = FALSE, CrI = CrI)
#
#            ## FOR LOOP HERE with number of correlations (corrs)
#            for(i in 1:corrs ) {
#                df <- array( NA, dim = c(standat$T, 3) )
#                df <- data.frame( rbind(df, forecasted_R[[i]]) )
#                df$period <- seq_len( dim(df )[1] ) 
#                df$Type <- factor( c( rep("Estimated", standat$T),
#                                     rep("Forecasted",  ahead)) )
#                
#                plt[[i]] <- retro$retro_plot[[i]] +
#                    geom_line(data = df, aes_string(x = "period",
#                                                    y = names(df)[1],
#                                                    color =  "Type") ) +                
#                    geom_errorbar(data = df,  aes_string(x = "period",
#                                                         ymin = names(df)[2],
#                                                         ymax = names(df)[3]),
#                                  inherit.aes = FALSE , alpha =  .3)+
#                    geom_point(data = df[df$Type == "Forecasted", ],
#                               mapping = aes_string(x = "period",
#                                                    y = names(df)[1],
#                                                    color =  "Type"), inherit.aes = FALSE) +
#                    coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
#                
#                
#                plot(plt[[i]])
#                if ( i == 1 ) {
#                    devAskNewPage( ask = TRUE)
#                }
#            }
#        }
#        ## print
#        return( forecasted_R )
#    }
#}

# TODO
print.forecast.bmgarch <- function(x, ...) {
    ahead <- x$forecast$meta$TS_length
    nt <- x$meta$nt
    TS_names <- x$meta$TS_names
    digits <- x$meta$digits

    # Mean structure
    if(x$meta$meanstructure == 1) {
        .sep()
        cat("[Mean]", "Forecast for", ahead, "ahead:")
        .newline(2)
        for(t in 1:nt) {
            cat(TS_names[t], ":")
            .newline()
            print(round(x$forecast$mean[,,t], digits))
        }
    }

    # Variance
    .sep()
    cat("[Variance]", "Forecast for", ahead, "ahead:")
    .newline(2)
    for(t in 1:nt) {
        cat(TS_names[t], ":")
        print(round(x$forecast$var[,,t], digits))
        .newline()
    }
    # Cors
    if(x$meta$param != "CCC") {
        cat("[Correlation]", "Forecast for", ahead, "ahead:")
        .newline(2)
        for(t in 1:nt) {
            cat(dimnames(x$forecast$cor)[[3]][t], ":")
            print(round(x$forecast$cor[,,t], digits))
            .newline()
        }
    }

    return(invisible(x))
}

# TODO Not sure if necessary.
summary.forecast.bmgarch <- function(object, ...) {
    
}

# TODO
plot.forecast.bmgarch <- function(x, ...) {
    
}

# Should create something tidy-like. Mainly, rbinding back/forecast, and adding columns for param, period, etc. 
as.data.frame.forecast.bmgarch <- function(x, ...) {
    
}

# Not sure if necessary (or feasible, since doing so requires mixed types.)
as.matrix.forecast.bmgarch <- function(x, ...) {
    
}

