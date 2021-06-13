##' @title Refit model
##' @param object bmgarch model object
##' @param data new data
##' @param xC_data new predictor
##' @keywords internal
.refit <- function(object,  data, xC_data  ) {
    fit_past <- bmgarch(data,
                        xC = xC_data,
                        parameterization = object$param,
                        P = object$mgarchP,
                        Q = object$mgarchQ,
                        chains = object$chains,
                        iterations = object$iter,
                        standardize_data = FALSE,
                        distribution = object$distribution,
                        meanstructure = object$meanstructure )
}

##' @keywords internal
## compute log of raw importance ratios
## sums over observations *not* over posterior samples
.sum_log_ratios <- function(ll, ids = NULL) {
if (!is.null(ids)) ll <- ll[, ids , drop = FALSE]
   - rowSums(ll)
}
 
##' @keywords internal
## more stable than log(sum(exp(x)))
.log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
 }
 
##' @keywords internal
.log_lik <- function(x) {
    rstan::extract(x$model_fit, pars = "log_lik")$log_lik
}

##' @keywords internal
## more stable than log(mean(exp(x))) 
.log_mean_exp <- function(x) {
   .log_sum_exp(x) - log(length(x))
}

##' \code{lfocv} returns the LFO-CV ELPD by either computing the exact ELDP or
##' by approximating it via
##' forward or backward approximation strategies based on Pareto smoothed
##' importance sampling
##' described in \insertCite{Buerkner2019}{bmgarch}.
##' @title Leave-Future-Out Cross Validation (LFO-CV)
##' @param x Fitted bmgarch model. \code{lfocv} inherits all attributes
##' from the bmgarch object
##' @param type Takes \code{lfo} (default) or \code{loo}. LFO-CV is recommended
##' for time-series but LOO-CV may be obtained to assess the structural part of the model.  
##' @param L Minimal length of times series before computing LFO
##' @param M M step head predictions. Defines to what period the LFO-CV should be tuned to. Defaults to M=1.  
##' @param mode backward elpd_lfo approximation, or exact elpd-lfo; 
##' Takes 'backward', and 'exact'. 'exact' fits N-L models and may
##' take a \emph{very} long time to complete. \code{forward} works too but is not
##' complete yet.
##' @param ... Not used
##' @return Approximate LFO-CV value and log-likelihood values across (L+1):N
##' timepoints
##' @references
##' \insertAllCited{}
##' @aliases loo
##' @importFrom loo loo
##' @importFrom stats sd weights
##' @examples
##' \dontrun{
##' data(stocks)
##' # Fit a DCC model 
##' fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
##'                parameterization = "DCC", standardize_data = TRUE,
##'                iterations = 500)
##'
##' # Compute expected log-predictive density (elpd) using the backward mode
##' # L is the upper boundary of the time-series before we engage in LFO-CV
##' lfob <- loo(fit, mode = 'backward',  L = 50 )
##' print(lfob)
##' }
##' @export 
##' @export loo
loo.bmgarch <- function(x, ..., type = 'lfo', L = NULL, M = 1, mode = "backward") {
    object <- x
    ahead <- M
    N <- object$TS_length
    ## List for returned objects
    outl <- list( )
    
    ## "Classic" loo, based on backcasted log_lik
    if( type == 'loo' ) {
        if( is.null( L ) ) L <- 0
        ll <- rstan::extract(object$model_fit, pars =  "log_lik")$log_lik
        LL <- ll[,  (L + 1):N]
        ## obtain chain id vector for relative_eff
        n_chains <- object$model_fit@sim$chains
        n_samples <- object$model_fit@sim$iter - object$model_fit@sim$warmup
        chain_id <- rep(seq_len(n_chains ),  each = n_samples )
        r_eff <- loo::relative_eff(exp(LL), chain_id = chain_id)
        backcast_loo <- loo::loo( LL, r_eff =  r_eff )
        outl$backcast_loo <- backcast_loo$estimates[,'Estimate']['elpd_loo']
        outl$type <- type
    } else if ( type == 'lfo' ) {
        if( is.null( L ) ) stop( "Provide L, length of time-series before fitting computing LFO-CV")
        
        ## Backward approach: start with log_lik from fully fitted model with all observations
        ## Go backwards until k threshold is reched and refit
        
        loglik <- matrix(nrow = dim( .log_lik( object ) )[1], ncol = N)
        approx_elpds_1sap <- rep(NA, N)
        ks <- NULL

        refits <- ks <- NULL
        
        if( mode == "backward" ) {
            k_thres <- 0.6 # Threshold of .6 based on Buerkner at al (2020) paper
            fit_past <- object
            i_refit <- N

            for (i in seq((N - ahead), L, by = -ahead) ) {
                loglik[, (i + 1):( i + ahead)] <- .log_lik(fit_past)[, (i + 1):( i + ahead)]
                logratio <- .sum_log_ratios(loglik, (i + 1):i_refit)
                psis_obj <- suppressWarnings(loo::psis(logratio))
                k <- loo::pareto_k_values(psis_obj)
                ks <- c(ks, k)
                if (k > k_thres) {
                    ## refit the model based on the first i observations
                    i_refit <- i
                    refits <- c(refits, i)
                    past <- 1:i
                    oos <- (i+1):(i+ahead)
                    df_past <- object$RTS_full[past, , drop = FALSE]
                    xC_past <- object$xC[past, , drop = FALSE]
                    df_oos <- object$RTS_full[c(past, oos), , drop = FALSE]
                    xC_oos <- object$xC[c(past, oos), , drop = FALSE]
                    fit_past <- .refit(object, data = df_past, xC_data = xC_past )
                    fc <- bmgarch::forecast(object = fit_past, ahead = ahead,
                                            xC = xC_oos[oos, ,drop = FALSE],
                                            newdata = df_oos[oos,,drop = FALSE])
                    loglik[, (i+1):(i+ahead) ] <- fc$forecast$log_lik[[1]]
                    if(ahead == 1 ) {
                        approx_elpds_1sap[ i+1 ] <- .log_mean_exp(loglik[, i+1 ])
                    } else {
                        approx_elpds_1sap[(i+1):(i+ahead)] <-
                            apply(loglik[, (i+1):(i+ahead) ], MARGIN = 2, FUN = .log_mean_exp )    
                    }
                } else {
                    lw <- weights(psis_obj, normalize = TRUE)[, 1]
                    if(ahead == 1 ) {
                        approx_elpds_1sap[ i+1 ] <-  .log_sum_exp(lw + loglik[, i+1 ])    
                    } else {
                        approx_elpds_1sap[(i+1):(i+ahead)] <-
                            apply((lw + loglik[, (i+1):(i+ahead)]),2,.log_sum_exp )
                    }
                }
            }
            out <- approx_elpds_1sap
            
        } else if ( mode == "forward" ) {
            warn <- function() warning("'forward' method not fully implemented")
            if( mode == 'forward' ) warn( )
            
            k_thres <- 0.7

            out <- rep( NA, N )
            exact_elpds_1sap <- rep( NA, N )

            df_dat <- object$RTS_full
            xC_dat <- object$xC
            
            oos <- L + 1
            
            ## Refit the model using the first L observations
            df_start <- object$RTS_full[1:L, , drop = FALSE]
            xC_start <- object$xC[1:L, , drop = FALSE]
            fit_start <- .refit(object, data = df_start, xC_data = xC_start )

            ## Exact ELPD: Computed from log_lik of forecast function
            fc <- bmgarch::forecast(fit_start, ahead = ahead,
                                    xC = xC_dat[ oos, , drop = FALSE],
                                    newdata = df_dat[ oos, ,drop = FALSE])

            ## Exact log_lik
            loglik[, oos ] <- fc$forecast$log_lik[[1]]
            exact_elpds_1sap[ L ] <- .log_mean_exp( loglik[, oos ] )
            out[ L ] <- exact_elpds_1sap[L]
            ## Also write content to approx_elpds as the summary for
            ## "forward" is obtained over all approx_elpds's
            approx_elpds_1sap[ L ] <- exact_elpds_1sap[L]
            i_refit <- L + 1

            if( L < N - ahead ) {
                for (i in (L + 1):(N-ahead) ) {
                    
                    logratio <- .sum_log_ratios(loglik, i_refit:i )#(L+1):(i+1))
                    psis_obj <- suppressWarnings(loo::psis(logratio))
                    k <- loo::pareto_k_values(psis_obj)
                    ks <- c(ks, k)
                    print( ks )
                    if( k > k_thres ) {
                        refits <- c(refits, i)

                        df_start <- object$RTS_full[1:i, , drop = FALSE]
                        xC_start <- object$xC[1:i, , drop = FALSE]

                        fit_start <- .refit(object, data = df_start, xC_data = xC_start )
                        ahead <- 1
                        fc <- bmgarch::forecast(fit_start, ahead = ahead,
                                                xC = xC_dat[ i+1, , drop = FALSE],
                                                newdata = df_dat[ i+1, ,drop = FALSE])

                        ## Exact log_lik
                        loglik[, i+1 ] <- fc$forecast$log_lik[[1]]
                        
                        exact_elpds_1sap[ i ] <- .log_mean_exp( loglik[, i+1 ] )
                        
                        logratio <- .sum_log_ratios(loglik, i+1 )#(L+1):(i+1))
                        out[ i ] <- exact_elpds_1sap[ i ]
                        psis_obj <- suppressWarnings(loo::psis(logratio))
                        k <- loo::pareto_k_values(psis_obj)
                        i_refit <- i+1
                                        #k <- 0
                    } else {
                        lw <- weights(psis_obj, normalize = TRUE)[, 1]
                        ahead <- ( i+1 )-( i_refit-1 )
                        fc <- bmgarch::forecast(fit_start, ahead = ahead,
                                                xC = xC_dat[ ( i_refit ):(i+1), , drop = FALSE],
                                                newdata = df_dat[  ( i_refit ):(i+1), ,drop = FALSE])

                        loglik[, i+1 ] <- fc$forecast$log_lik[[1]][, ahead]
                        approx_elpds_1sap[ i ] <-  .log_sum_exp(lw + loglik[, i+1 ])
                        out[i] <- approx_elpds_1sap[ i ] 
                    }          
                }
            }
        } else if (mode == "exact" ) {
            k_thres <- 0
            refits <- N - ( L+1 )
            df_dat <- object$RTS_full
            xC_dat <- object$xC
            exact_elpds_1sap <- rep( NA, N )

            for(i in L:( N-1 ) ) {
                fit_start <- .refit(object, data = df_dat[1:i, , drop = FALSE], xC_data = xC_dat[1:i, , drop = FALSE] )
                fc <- bmgarch::forecast(fit_start, ahead = ahead, xC = xC_dat[ i+1, , drop = FALSE], newdata = df_dat[ i+1, ,drop = FALSE])
                loglik[, i+1 ] <-  fc$forecast$log_lik[[1]]
            }
            loglik_exact <- loglik[, (L+1):N]
            exact_elpds_1sap <- apply(loglik_exact, 2, .log_mean_exp )
            out <- exact_elpds_1sap
        } else {
            stop("'mode' needs to be either 'forward', 'backward' or 'exact'." )
        }

        
        refit_info <- cat("Using threshold ", k_thres, 
                          ", model was refit ", length(refits), 
                          " times, at observations", refits, "\n")
    
        ## Return relevant objects
        
        outl$refits <-refits
        outl$ks <- ks
        outl$approx_elpd_1sap <- sum( approx_elpds_1sap, na.rm = TRUE )
        outl$out <- out
        outl$mode <- mode
        outl$L <- L
        outl$type <- type
        
        outl$loglik <- loglik[, ( L + 1 ):N ]
    }
    
    attr(outl, "class" ) <- "loo.bmgarch"
    return( outl )
}

##' @title print method for lfocv
##' @param x lfo object
##' @param ... Not used.
##' @return Invisible lfocv object
##' @author philippe
##' @export
print.loo.bmgarch <- function( x, ... ) {
    if( x$type == 'loo' ) {
        cat('elpd_loo ', x$backcast_loo)
    } else if( x$type == 'lfo' ) {
        if(x$mode == 'backward' | x$mode == 'forward' ) {
            cat("Model was fit using ", x$mode, " mode.\n")
            cat("Approximate ELPD_LFO: ", x$approx_elpd_1sap, "\n" )
        } else if(x$mode == 'exact' ){
            cat("Model was fit using ", x$mode, " mode.\n")
            cat("Exact ELPD_LFO: ", sum( x$out, na.rm = TRUE), "\n" )
        }
    }
    return(invisible(x))
}
