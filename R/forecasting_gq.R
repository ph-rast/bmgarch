##' Estimates (weighted) forecasted means, variances, and correlations from a fitted bmgarch model.
##' @title Forecast method for bmgarch objects.
##' @param object bmgarch object.
##' @param ahead Integer (Default: 1). Periods to be forecasted ahead.
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c. Used in a log-linear model on the constant variance terms. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param newdata Future datapoints for LFO-CV computation
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param seed Integer (Optional). Specify seed for \code{\link[rstan]{sampling}}.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @param weights Takes weights from model_weight function. Defaults to 1 -- this parameter is not typically set by user.
##' @param L Minimal length of time series before engaging in lfocv 
##' @param method Ensemble methods, 'stacking' (default) or 'pseudobma'
##' @param inc_samples Logical (Default: FALSE). Whether to return the MCMC samples for the fitted values.
##' @param ... Not used 
##' @return forecast.bmgarch object. List containing \code{forecast}, \code{backcast}, and \code{meta}data.
##' See \code{\link{fitted.bmgarch}} for information on \code{backcast}.
##' \code{forecast} is a list of four components:
##' \describe{
##'   \item{mean}{\code{[N, 7, TS]} array of mean forecasts, where N is the timeseries length, and TS is the number of time series. E.g., \code{fc$forecast$mean[3,,"tsA"]} is the 3-ahead mean forecast for time series "tsA".}
##'   \item{var}{\code{[N, 7, TS]} array of variance forecasts, where N is the timeseries length, and TS is the number of time series. E.g., \code{fc$forecast$var[3,,"tsA"]} is the 3-ahead variance forecast for time series "tsA".}
##'   \item{cor}{\code{[N, 7, TS(TS - 1)/2]} array of correlation forecasts, where N is the timeseries length, and \code{TS(TS - 1)/2} is the number of correlations. E.g., \code{fc$forecast$cor[3,, "tsB_tsA"]} is the 3-ahead forecast for the correlation between "tsB" and "tsA". Lower triangular correlations are saved.}
##'   \item{meta}{Meta-data specific to the forecast. I.e., TS_length (number ahead) and xC.}
##'   \item{samples}{List}. If inc_samples is \code{TRUE}, then a list of arrays of MCMC samples for means, vars, and cors. Each array is [Iteration, Period, ..., ...].
##' }
##' @aliases forecast
##' @importFrom forecast forecast
##' @export 
##' @export forecast
##' @examples
##' \dontrun{
##' data(panas)
##' # Fit DCC(2,2) with constant mean structure.
##' fit <- bmgarch(panas, parameterization = "DCC", P = 2, Q = 2, meanstructure = "constant")
##'
##' # Forecast 8 ahead
##' fit.fc <- forecast(fit, ahead = 8)
##'
##' # Print forecasts
##' fit.fc
##' print(fit.fc)
##'
##' # Plot variance forecasts
##' plot(fit.fc, type = "var")
##'
##' # Plot correlation forecasts
##' plot(fit.fc, type = "cor")
##'
##' # Save backcasted and forecasted values as data frame.
##' fit.fc.df <- as.data.frame(fit.fc)
##'
##' # Save only forecasted values as data frame.
##' fit.fc.df <- as.data.frame(fit.fc, backcast = FALSE)
##'
##' # Add another model, compute model weights and perform a model weighted forecast
##'
##' # Fit a DCC(1,1) model
##' fit1 <- bmgarch(panas, parameterization = "DCC", P = 1, Q = 1, meanstructure = "constant")
##'
##' # Compute model stacking weights based on the last 19 time points (with L = 80)
##' blist <- bmgarch_list( fit1, fit )
##' mw <- model_weights(blist, L = 80)
##'
##' # Weighted forecasts:
##' w.fc <- forecast(object = blist, ahead = 8, weights = mw)
##' }
forecast.bmgarch <- function(object, ahead = 1, xC = NULL,
                             newdata = NULL, CrI = c(.025, .975),
                             seed = NA, digits = 2, weights = NULL,
                             L = NA, method = 'stacking', inc_samples = FALSE, ...) {
    
    ## Are we dealing with one object or a list of objects
    n_mods <- 1
    if("bmgarch_list" %in% class(object)) {
        n_mods <- length(object)
    } else {
        object <- bmgarch_list(object)
    }
    TS_names <- object[[1]]$TS_names
    # Check for TS name consistency
    TS_names_consistent <- all(sapply(object, function(x) {
                                all(x$TS_names == TS_names)
                            }))
    if(!TS_names_consistent) {
        # Could *possibly* rearrange the column orders to 'fix' this, but this is a much safer default.
        # Could also check whether the training data are the same across models, to ensure the predictions make sense.
        stop("Time series column names are not consistent across models. Forecasting halted.")
    }
    ## n_mods <- .depth( object )
    # Define a 0 array for stan.
    if(is.null(xC)) {
        xC <- array(0, dim = c(ahead, object[[1]]$nt))
    }
    if(is.null(newdata)) {
        newdata <- array(0, dim = c(ahead, object[[1]]$nt))
        compute_log_lik <- 0
    } else {
        compute_log_lik <- 1
    }

    ## if user provides weights from the model_weigths function
    ## proceed directly to forecasting, else, run model_weights
    ## and extract model weights here
    ## Case 1: No model weights provided
    ## Case 2: model weights from a model_weight object
    ## Case 3: No model weights requested
    if(n_mods > 1 & is.null( weights ) ) {
        mw <- bmgarch::model_weights(bmgarch_objects = object, L = L)
        weights <- mw$wts[]
    } else if( n_mods > 1 & !is.null( weights ) ) {
        weights <- weights$wts[]
    } else if( n_mods == 1 ) {
        weights <- 1
        ## object[[1]] <- object
    }
    
    object.f <- lapply(object, function(m) {
        standat <- list(T = m$TS_length,
                        nt = m$nt,
                        rts = cbind(m$RTS_full),
                        xC = m$xC,
                        Q =  m$mgarchQ,
                        P =  m$mgarchP,
                        ahead =  ahead, 
                        meanstructure =  m$meanstructure,
                        distribution =  m$num_dist,
                        xC_p =  xC,
                        future_rts = newdata,
                        compute_log_lik =  compute_log_lik)

        gqs_model <- switch(m$param,
                            DCC = stanmodels$forecastDCC,
                            CCC = stanmodels$forecastCCC,
                            BEKK =stanmodels$forecastBEKK,
                            pdBEKK = stanmodels$forecastBEKK,
                            NULL)
        if(is.null(gqs_model)) {
            stop("bmgarch object 'param' does not match a supported model. ",
                    m$param, "is not one in ", paste0(supported_models, collapse = ", "),
                    ".")
        }
        backcast <- max(m$mgarchP, m$mgarchQ)
        nt <- m$nt
        cast_start <- (m$TS_length - backcast + 1)
        forecast_start <- (m$TS_length + 1)
        forecast_end <- (m$TS_length + ahead)

        ## TODO: Limit pars to only what is needed (H_p, R/R_p, rts_p, mu_p)
        forecasted <- rstan::gqs(gqs_model,
                                 draws = as.matrix(m$model_fit),
                                 data = standat,
                                 seed = seed)
        return(forecasted)
    })

    ## Init f.mean
    f.mean <- .get_stan_summary(object.f, "rts_forecasted", CrI, weights)

    ## f.var
    f.var <- .get_stan_summary(object.f, "H_forecasted", CrI, weights)

    ## Init f.cor
    f.cor <- .get_stan_summary(object.f, "R_forecasted", CrI, weights)

    # Restructure to array
    ## backcast <- max(object[[1]]$mgarchP, object[[1]]$mgarchQ)
    nt <- object[[1]]$nt
    ## cast_start <- (object[[1]]$TS_length - backcast + 1)
    forecast_start <- (object[[1]]$TS_length + 1)
    forecast_end <- (object[[1]]$TS_length + ahead)

    ## f.mean
    stan_sum_cols <- colnames(f.mean)
    f.mean <- array(f.mean, dim = c(nt, ahead, ncol(f.mean)))
    f.mean <- aperm(f.mean, c(2,3,1))
    dimnames(f.mean) <- list(period = forecast_start:forecast_end, stan_sum_cols, TS = TS_names)

    ## f.var
    ### Pull out indices for [period, a, a]
    f.var.indices <- grep("H_forecasted\\[[[:digit:]]+,([[:digit:]]+),\\1]", rownames(f.var), value = TRUE)
    f.var <- f.var[f.var.indices,]
    f.var <- array(f.var, dim = c(nt, ahead, ncol(f.var)))
    f.var <- aperm(f.var, c(2, 3, 1))
    dimnames(f.var) <- list(period = forecast_start:forecast_end, stan_sum_cols, TS = TS_names)

    ## f.cor
    # Lower-triangular indices
    f.cor.indices.L <- which(lower.tri(matrix(0, nt, nt)), arr.ind = TRUE)
    # Labels mapping to TS names
    f.cor.indices.L.labels <- paste0(TS_names[f.cor.indices.L[, 1]], "_",
                                     TS_names[f.cor.indices.L[, 2]])
    # Indices as "a,b"
    f.cor.indices.L.char <- paste0(f.cor.indices.L[, 1], ",", f.cor.indices.L[,2])
    # Indicices as "[period,a,b]"
    f.cor.indices.L.all <- paste0("R_forecasted[",1:(ahead), ",",
                                  rep(f.cor.indices.L.char, each = (ahead)),"]")
    # Get only these elements.
    f.cor <- f.cor[f.cor.indices.L.all, ,drop = FALSE]
    f.cor <- array(f.cor, dim = c(ahead, length(f.cor.indices.L.char), ncol(f.cor ) ))
    f.cor <- aperm(f.cor, c(1, 3, 2))
    dimnames(f.cor) <-  list(period = forecast_start:forecast_end, stan_sum_cols, TS = f.cor.indices.L.labels)

    # Remove backcasts from forecasts.
    ## f.mean <- f.mean[-c(1:backcast), , , drop = FALSE]
    ## f.var <- f.var[-c(1:backcast), , , drop = FALSE]
    ## f.cor <- f.cor[-c(1:backcast), , , drop = FALSE]
   
    out <- list()
    out$forecast$mean <- f.mean
    out$forecast$var <- f.var
    out$forecast$cor <- f.cor
    out$forecast$meta <- list(xC = xC, TS_length = ahead)

    if(inc_samples) {
        out$forecast$samples$mean <- .weighted_samples(object.f, "rts_forecasted", weights)$rts_forecasted[, ,, drop = FALSE]
        out$forecast$samples$var <- .weighted_samples(object.f, "H_forecasted", weights)$H_forecasted[,, , , drop = FALSE]
        out$forecast$samples$cor <- .weighted_samples(object.f, "R_forecasted", weights)$R_forecasted[,, , , drop = FALSE]
        ## out$forecast$samples$mean <- out$forecast$samples$mean[] # Todo, remove backcast
    }

        ## Extract all log_lik simulations
    if(compute_log_lik == 1 ) {
        log_lik <- lapply(object.f, function(x) {
            rstan::extract(x, pars = "log_lik")$log_lik
        })
        out$forecast$log_lik <- log_lik
    }

    
    metaNames <- c("param", "distribution", "num_dist", "nt", "TS_length", "TS_names", "RTS_full", "mgarchQ", "mgarchP", "xC", "meanstructure")
    meta <- with(object[[1]], mget(metaNames))
    meta_bmgarch_list <- lapply(object, function(x) {with(x, mget(metaNames))})
    out$meta_list <- meta_bmgarch_list
    out$meta <- meta
    out$meta$n_mods <- n_mods
    out$meta$digits <- digits
    out$meta$CrI <- CrI
    out$meta$weights <- weights

    out$backcast <- fitted.bmgarch(object, CrI, digits = digits, weights = weights, inc_samples = inc_samples)$backcast

    class(out) <- "forecast.bmgarch"
    return(out)
}
##' Extracts the model-predicted means, variances, and correlations for the fitted data.
##'
##' Whereas \code{\link{forecast.bmgarch}} computes the \emph{forecasted} values for future time periods, \code{fitted.bmgarch} computes the \emph{backcasted} (model-predicted) values for the observed time periods.
##' @title Fitted (backcasting) method for bmgarch objects.
##' @param object bmgarch object.
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @param weights Takes weights from model_weight function. Defaults to 1 -- this parameter is not typically set by user.
##' @param inc_samples Logical (Default: FALSE). Whether to return the MCMC samples for the fitted values.
##' @param ... Not used.
##' @return fitted.bmgarch object. List containing \code{meta}data and the \code{backcast}. Backcast is a list containing three elements:
##' \describe{
##'   \item{mean}{\code{[N, 7, TS]} array of mean backcasts, where N is the timeseries length, and TS is the number of time series. E.g., \code{bc$backcast$mean[3,,"tsA"]} is the mean backcast for the third observation in time series "tsA".}
##'   \item{var}{\code{[N, 7, TS]} array of variance backcasts, where N is the timeseries length, and TS is the number of time series. E.g., \code{bc$backcast$var[3,,"tsA"]} is the variance backcast for the third observation in time series "tsA".}
##'   \item{cor}{\code{[N, 7, TS(TS - 1)/2]} array of correlation backcasts, where N is the timeseries length, and \code{TS(TS - 1)/2} is the number of correlations. E.g., \code{bc$backcast$cor[3,, "tsB_tsA"]} is the backcast for the correlation between "tsB" and "tsA" on the third observation. Lower triangular correlations are saved.}
##'   \item{samples}{List}. If inc_samples is \code{TRUE}, then a list of arrays of MCMC samples for means, vars, and cors. Each array is [Iteration, Period, ..., ...].
##' }
##' @importFrom stats fitted
##' @export
##' @examples
##' \dontrun{
##' data(panas)
##' # Fit CCC(1,1) and constant meanstructure.
##' fit <- bmgarch(panas, parameterization = "CCC", meanstructure = "constant")
##'
##' # Obtain fitted values
##' fit.bc <- fitted(fit)
##'
##' # Print fitted values
##' print(fit.bc)
##'
##' # Plot fitted values (plot.bmgarch calls fitted internally)
##' plot(fit, type = "var")
##'
##' # Save fitted values as data frame
##' fit.bc.df <- as.data.frame(fit.bc)
##' }
fitted.bmgarch <- function(object, CrI = c(.025, .975), digits = 2, weights = NULL, inc_samples = FALSE, ...) {
    n_mods <- 1
    if("bmgarch_list" %in% class(object)) {
        n_mods <- length(object)
    } else {
        object <- bmgarch_list(object)
    }
    nt <- object[[1]]$nt
    TS_length <- object[[1]]$TS_length
    TS_names <- object[[1]]$TS_names

    if(n_mods > 1 & is.null(weights)) {
        stop("Weights must be provided.")
    } else if(n_mods == 1) {
        weights <- 1
    }

    fits <- lapply(object, function(m) {m$model_fit})
    b.mean <- .get_stan_summary(fits, "mu", CrI, weights)
    b.var <- .get_stan_summary(fits, "H", CrI, weights)
    b.cor <- .get_stan_summary(fits, "corH", CrI, weights)

    # Restructure
    ## b.mean
    stan_sum_cols <- colnames(b.mean)
    b.mean <- array(b.mean, dim = c(nt, TS_length, ncol(b.mean)))
    b.mean <- aperm(b.mean, c(2,3,1))
    dimnames(b.mean) <- list(period = 1:TS_length, stan_sum_cols, TS = TS_names)

    ## b.var
    b.var.indices <- grep("H\\[[[:digit:]]+,([[:digit:]]+),\\1]", rownames(b.var), value = TRUE)
    b.var <- b.var[b.var.indices,]
    b.var <- array(b.var, dim = c(nt, TS_length, ncol(b.var)))
    b.var <- aperm(b.var, c(2, 3, 1))
    dimnames(b.var) <- list(period = 1:TS_length, stan_sum_cols, TS = TS_names)

    ## b.cor
    # Lower-triangular indices
    b.cor.indices.L <- which(lower.tri(matrix(0, nt, nt)), arr.ind = TRUE)
    # Labels mapping to TS names
    b.cor.indices.L.labels <- paste0(TS_names[b.cor.indices.L[,1]], "_", TS_names[b.cor.indices.L[,2]])
    # Indices as "a,b"
    b.cor.indices.L.char <- paste0(b.cor.indices.L[,1], ",", b.cor.indices.L[,2])
    # Indicices as "[period,a,b]"
    b.cor.indices.L.all <- paste0("corH[",1:TS_length, ",", rep(b.cor.indices.L.char, each = TS_length),"]")
    # Get only these elements.
    b.cor <- b.cor[b.cor.indices.L.all,]
    b.cor <- array(b.cor, dim = c(TS_length, length(b.cor.indices.L.char), ncol(b.cor)))
    b.cor <- aperm(b.cor, c(1, 3, 2))
    dimnames(b.cor) <- list(period = 1:TS_length, stan_sum_cols, TS = b.cor.indices.L.labels)

    out <- list()
    out$backcast$mean <- b.mean
    out$backcast$var <- b.var
    out$backcast$cor <- b.cor

    if(inc_samples) {
        out$backcast$samples$mean <- .weighted_samples(fits, "mu", weights)$mu
        out$backcast$samples$var <- .weighted_samples(fits, "H", weights)$H
        out$backcast$samples$cor <- .weighted_samples(fits, "corH", weights)$corH
    }

    metaNames <- c("param", "distribution", "num_dist", "nt", "TS_length", "TS_names", "RTS_full", "mgarchQ", "mgarchP", "xC", "meanstructure")
    meta <- with(object[[1]], mget(metaNames))
    out$meta_list <- lapply(object, function(x) {with(x, mget(metaNames))})
    out$meta <- meta
    out$meta$digits <- digits
    out$meta$n_mods <- n_mods
    out$meta$CrI <- CrI
    out$meta$weights <- weights

    class(out) <- "fitted.bmgarch"
    return(out)
}


##' @title Print method for forecast.bmgarch objects.
##' @param x forecast.bmgarch object. See \code{\link{forecast.bmgarch}}
##' @param ... Not used.
##' @return x (invisible).
##' @author Stephen R. Martin
##' @export
print.forecast.bmgarch <- function(x, ...) {
    ahead <- x$forecast$meta$TS_length
    nt <- x$meta$nt
    TS_names <- x$meta$TS_names
    digits <- x$meta$digits

    if(x$meta$n_mods > 1) {
        .sep()
        cat("LFO-weighted forecasts across ", x$meta$n_mods, "models.")
        .newline()
    }

    # Mean structure
    meanstructure <- any(sapply(x$meta_list, function(x) {x$meanstructure == 1}))
    ## if(x$meta$meanstructure == 1 | x$meta$n_mod > 1) {
    if(meanstructure) {
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
        .newline()
        print(round(x$forecast$var[,,t], digits))
    }
    # Cors
    condCor <- any(sapply(x$meta_list, function(x) {x$param != "CCC"}))
    ## if(x$meta$param != "CCC" | x$meta$n_mod > 1) {
    if(condCor) {
        cat("[Correlation]", "Forecast for", ahead, "ahead:")
        .newline(2)
        for(t in 1:(nt*(nt - 1) / 2)) {
            cat(dimnames(x$forecast$cor)[[3]][t], ":")
            .newline()
            print(round(x$forecast$cor[,,t], digits))
        }
    }


}

##' @title Print method for fitted.bmgarch objects.
##' @param x fitted.bmgarch object.
##' @param ... Not used.
##' @return object (invisible).
##' @author Stephen R. Martin
##' @export
print.fitted.bmgarch <- function(x, ...) {
    TS_length <- x$meta$TS_length
    nt <- x$meta$nt
    TS_names <- x$meta$TS_names
    digits <- x$meta$digits

    # Mean structure
    if(x$meta$meanstructure == 1) {
        .sep()
        cat("[Mean]", "Fitted values:")
        .newline(2)
        for(t in 1:nt) {
            cat(TS_names[t], ":")
            .newline()
            print(round(x$backcast$mean[,,t], digits))
        }
    }

    # Variance
    .sep()
    cat("[Variance]", "Fitted values:")
    .newline(2)
    for(t in 1:nt) {
        cat(TS_names[t], ":")
        .newline()
        print(round(x$backcast$var[,,t], digits))
    }

    # Cors
    if(x$meta$param != "CCC") {
        cat("[Correlation]", "Fitted values:")
        .newline(2)
        for(t in 1:(nt*(nt - 1) / 2)) {
            cat(dimnames(x$backcast$cor)[[3]][t], ":")
            .newline()
            print(round(x$backcast$cor[,,t], digits))
        }
    }
    object <- x
    return(invisible(object))
}

##' Helper function for as.data.frame.{fitted, forecast}. Converts predictive array to data.frame.
##' 
##' 
##' @title Convert predictive array to data.frame.
##' @param arr Array to convert into data frame.
##' @param type String. "backcast" or "forecast".
##' @param param String. "var", "mean", or "cor".
##' @return data.frame. Columns: period, type (backcast, forecast), param (var, mean, cor), TS (which time series, or which correlation for param = cor), summary columns.
##' @author Stephen R. Martin
##' @keywords internal
.pred_array_to_df <- function(arr, type = "backcast", param = "var") {
    dims <- dim(arr)
    arrnames <- dimnames(arr)

    dfList <- apply(arr, 3, function(x) {
        out <- as.data.frame(x)
        out$period <- as.numeric(rownames(x))
        out
    })
    for(i in seq_len(length(dfList))) {
        dfList[[i]]$TS <- arrnames[[3]][i]
    }
    df <- do.call(rbind, dfList)
    df$type <- type
    df$param <- param

    rownames(df) <- NULL

    return(df)
}

##' @title as.data.frame method for forecast.bmgarch objects.
##' @param x forecast.bmgarch object.
##' @param backcast Logical (Default: True). Whether to include "backcasted" values from \code{\link{fitted.bmgarch}} in data frame.
##' @param ... Not used.
##' @return Data frame.
##' @author Stephen R. Martin
##' @export
as.data.frame.forecast.bmgarch <- function(x, ..., backcast = TRUE ) {

    # Forecast
    dfList <- list()
    dfList$forecast.mean <- .pred_array_to_df(x$forecast$mean, "forecast", "mean")

    dfList$forecast.var <- .pred_array_to_df(x$forecast$var, "forecast", "var")

    ## if(x$meta$param != "CCC") {
    condCor <- any(sapply(x$meta_list, function(x) {x$param != "CCC"}))
    if(condCor) {
        dfList$forecast.cor <- .pred_array_to_df(x$forecast$cor, "forecast", "cor")
    }

    if(backcast) {
        # Backcast
        dfList$backcast.mean <- .pred_array_to_df(x$backcast$mean, "backcast", "mean")

        dfList$backcast.var <- .pred_array_to_df(x$backcast$var, "backcast", "var")

        if(x$meta$param != "CCC") {
            dfList$backcast.cor <- .pred_array_to_df(x$backcast$cor, "backcast", "cor")
        }
    }

    # Combine
    df <- do.call(rbind, dfList)

    # Re-order columns: period TS | type | param
    desc <- c("period","TS","type","param")
    cn <- colnames(df)
    cn_not_desc <- cn[!(cn %in% desc)]
    df <- df[,c(desc, cn_not_desc)]

    # Sort
    df <- df[with(df, order(param, TS, period)),]

    rownames(df) <- NULL

    return(df)

}

##' @title as.data.frame method for fitted.bmgarch objects.
##' @param x fitted.bmgarch object.
##' @param ... Not used.
##' @return Data frame.
##' @author Stephen R. Martin
##' @export
as.data.frame.fitted.bmgarch <- function(x, ...) {
    dfList <- list()

    dfList$backcast.mean <- .pred_array_to_df(x$backcast$mean, "backcast", "mean")

    dfList$backcast.var <- .pred_array_to_df(x$backcast$var, "backcast", "var")

    if(x$meta$param != "CCC") {
        dfList$backcast.cor <- .pred_array_to_df(x$backcast$cor, "backcast", "cor")
    }

    # Combine
    df <- do.call(rbind, dfList)

    # Re-order columns: period TS | type | param
    desc <- c("period","TS","type","param")
    cn <- colnames(df)
    cn_not_desc <- cn[!(cn %in% desc)]
    df <- df[,c(desc, cn_not_desc)]

    # Sort
    df <- df[with(df, order(param, TS, period)),]

    rownames(df) <- NULL

    return(df)

}

##' @title Collect bmgarch objects into list.
##' @param ... bmgarch objects.
##' @return List of bmgarch objects. Class: bmgarch_list and bmgarch.
##' @export
bmgarch_list <- function(...) {
    out <- list(...)
    class(out) <- c("bmgarch_list", "bmgarch")
    return(out)
}
