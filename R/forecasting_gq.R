##' Figure out whther object is a list of models or just one model
##' Check for nesting structure; If depth == 1, one object else list of models
##' from: \url{https://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r}
##' @title Obtain nesting depth of list
##' @param this 
##' @param thisdepth 
##' @return Depth
##' @author philippe
##' @keywords internal
.depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,.depth,thisdepth=thisdepth+1))))    
  }
}
##' Obtain SD's over columns in lists
##' @title Column SD's
##' @param x 
##' @return SD's at the columns level witin list
##' @author philippe
##' @keywords internal 
.colSDs <- function(x) {
    lapply(x, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), sd)
    })
}

##' Obtain quantiles over columns in lists
##' @title Quantiles within lists
##' @param x 
##' @param probs Quantile(s). Inherits from \code{forecast} which defaults to \code{c(.025, .975)}.
##' @return Quantiles at the column level within lists
##' @author philippe
##' @keywords internal
.colQTs <- function(x, probs = CrI) {
    lapply(x, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), quantile, probs)
    })
}

##' @title Sort list object 
##' @param x 
##' @return sorted list as vector
##' @author philippe
##' @keywords internal
.sort <- function(x ) {
    c(apply( x, 1, FUN = function(x) {
        c(x )
    }))
}



##' Estimates forecasted means, variances, and correlations from a fitted bmgarch model.
##' @title Forecast method for bmgarch objects.
##' @param object bmgarch object.
##' @param ahead Integer (Default: 1). Periods to be forecasted ahead.
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c. Used in a log-linear model on the constant variance terms. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param newdata Future datapoints for LFO-CV computation
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param seed Integer (Optional). Specify seed for \code{\link[rstan]{sampling}}.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @param weights Takes weights from model_weight function. Defaults to 1 -- this parameter is not typically set by user.
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
##' @author Stephen R. Martin
##' @importFrom forecast forecast
##' @export forecast
##' @export
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
##' }
forecast.bmgarch <- function(object, ahead = 1, xC = NULL,
                             newdata = NULL, CrI = c(.025, .975),
                             seed = NA, digits = 2, weights = NULL,
                             L = NA, method = 'stacking', inc_samples = FALSE) {
    
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
                            DCC = bmgarch:::stanmodels$forecastDCC,
                            CCC = bmgarch:::stanmodels$forecastCCC,
                            BEKK = bmgarch:::stanmodels$forecastBEKK,
                            pdBEKK = bmgarch:::stanmodels$forecastBEKK,
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
    f.cor <- f.cor[f.cor.indices.L.all,]
    f.cor <- array(f.cor, dim = c(ahead, length(f.cor.indices.L.char), 7 ))
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
##' @author Stephen R. Martin
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
##' @param object fitted.bmgarch object.
##' @param ... Not used.
##' @return object (invisible).
##' @author Stephen R. Martin
##' @export
print.fitted.bmgarch <- function(object, ...) {
    TS_length <- object$meta$TS_length
    nt <- object$meta$nt
    TS_names <- object$meta$TS_names
    digits <- object$meta$digits

    # Mean structure
    if(object$meta$meanstructure == 1) {
        .sep()
        cat("[Mean]", "Fitted values:")
        .newline(2)
        for(t in 1:nt) {
            cat(TS_names[t], ":")
            .newline()
            print(round(object$backcast$mean[,,t], digits))
        }
    }

    # Variance
    .sep()
    cat("[Variance]", "Fitted values:")
    .newline(2)
    for(t in 1:nt) {
        cat(TS_names[t], ":")
        .newline()
        print(round(object$backcast$var[,,t], digits))
    }

    # Cors
    if(object$meta$param != "CCC") {
        cat("[Correlation]", "Fitted values:")
        .newline(2)
        for(t in 1:(nt*(nt - 1) / 2)) {
            cat(dimnames(object$backcast$cor)[[3]][t], ":")
            .newline()
            print(round(object$backcast$cor[,,t], digits))
        }
    }

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
as.data.frame.forecast.bmgarch <- function(x, backcast = TRUE, ...) {

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
##' @author Stephen R. Martin
##' @export
bmgarch_list <- function(...) {
    out <- list(...)
    class(out) <- c("bmgarch_list", "bmgarch")
    return(out)
}
