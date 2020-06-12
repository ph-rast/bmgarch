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

##' \code{lfocv} returns the LFO-CV ELPD by either computing the exact ELDP or by approximating it via
##' forward or backward approximation strategies based on Pareto smoothed importance sampling
##' described in \insertCite{Buerkner2019}{bmgarch}. 
##' @title Leave-Future-Out Cross Validation (LFO-CV)
##' @param object Fitted bmgarch model. \code{lfocv} inherits all attributes from the bmgarch object.
##' @param L Minimal length of times series before computing LFO 
##' @param mode forward or backward elpd_lfo approximation, or exact elpd-lfo.
##' Takes 'forward', 'backward', and 'exact'. 'exact' fits N-L models and may
##' take a \emph{very} long time to complete. 
##' @return Approximate LFO-CV value and log-likelihood values across (L+1):N timepoints
##' @author philippe
##' @importFrom Rdpack reprompt
##' @references
##' \insertAllCited{}
##' 
##' @export
lfocv <- function(object, L = 50, mode = "backward") {
    m <- 1 ## currently on only lfo-1-ahead
    ahead <- m
    N <- object$TS_length

    ## Backward approach: start with log_lik from fully fitted model with all observations
    ## Go backwards until k threshold is reched and refit
    
    loglik <- matrix(nrow = dim( .log_lik( object ) )[1], ncol = N)
    approx_elpds_1sap <- rep(NA, N)
    ks <- NULL

    refits <- ks <- NULL
    i_refit <- N
    
    if( mode == "backward" ) {
        k_thres <- 0.6
        fit_past <- object
        

        for (i in (N - 1):L) {
            loglik[, i + 1] <- .log_lik(fit_past)[, i + 1]
            logratio <- .sum_log_ratios(loglik, (i + 1):i_refit)
            psis_obj <- suppressWarnings(loo::psis(logratio))
            k <- loo::pareto_k_values(psis_obj)
            ks <- c(ks, k)
            if (k > k_thres) {
                ## refit the model based on the first i observations
                i_refit <- i
                refits <- c(refits, i)
                past <- 1:i
                oos <- i + ahead ## come back and fix for when m > 1
                df_past <- object$RTS_full[past, , drop = FALSE]
                xC_past <- object$xC[past, , drop = FALSE]
                df_oos <- object$RTS_full[c(past, oos), , drop = FALSE]
                xC_oos <-  object$xC[c(past, oos), , drop = FALSE]
                fit_past <- .refit(object, data = df_past, xC_data = xC_past )
                fc <- bmgarch::forecast(fit_past, ahead = ahead, xC = xC_oos[oos, ,drop = FALSE],
                                        newdata = df_oos[oos,,drop = FALSE])
                loglik[, i + 1 ] <- fc$forecast$log_lik
                approx_elpds_1sap[i + 1] <- .log_mean_exp(loglik[, i + 1])
            
            } else {
                lw <- weights(psis_obj, normalize = TRUE)[, 1]
                approx_elpds_1sap[i + 1] <-  .log_sum_exp(lw + loglik[, i + 1])
            }
        }
        out <- approx_elpds_1sap
  
    } else if ( mode == "forward" ) {
        k_thres <- 0.7

        out <- rep( NA, N )
        exact_elpds_1sap <- rep( NA, N )

        df_dat <- object$RTS_full
        xC_dat <- object$xC
        
        oos <- L + ahead
        
        ## Refit the model using the first L observations
        df_start <- object$RTS_full[1:L, , drop = FALSE]
        xC_start <- object$xC[1:L, , drop = FALSE]
        fit_start <- .refit(object, data = df_start, xC_data = xC_start )

        oos
        ## Exact ELPD: Computed from log_lik of forecast function
        fc <- bmgarch::forecast(fit_start, ahead = ahead, xC = xC_dat[ oos, , drop = FALSE], newdata = df_dat[ oos, ,drop = FALSE])

        ## Exact log_lik
        loglik[, oos ] <- fc$forecast$log_lik
        exact_elpds_1sap[ L ] <- .log_mean_exp( loglik[, oos ] )
        out[ L ] <- exact_elpds_1sap[L]    

        for (i in (L + 1):(N-ahead) ) {
            logratio <- .sum_log_ratios(loglik, i )#(L+1):(i+1))
            psis_obj <- suppressWarnings(loo::psis(logratio))
            k <- loo::pareto_k_values(psis_obj)
            ks <- c(ks, k)
            print( ks )
            if( k > k_thres ) {
                i_refit <- i
                refits <- c(refits, i)

                df_start <- object$RTS_full[1:i, , drop = FALSE]
                xC_start <- object$xC[1:i, , drop = FALSE]

                fit_start <- .refit(object, data = df_start, xC_data = xC_start )

                fc <- bmgarch::forecast(fit_start, ahead = ahead, xC = xC_dat[ i+1, , drop = FALSE], newdata = df_dat[ i+1, ,drop = FALSE])

                ## Exact log_lik
                loglik[, i+1 ] <- fc$forecast$log_lik
                
                exact_elpds_1sap[ i ] <- .log_mean_exp( loglik[, i+1 ] )
                
                logratio <- .sum_log_ratios(loglik, i+1 )#(L+1):(i+1))
                out[ i ] <- exact_elpds_1sap[ i ]
                psis_obj <- suppressWarnings(loo::psis(logratio))
                k <- loo::pareto_k_values(psis_obj)
                #k <- 0
            } else {
                lw <- weights(psis_obj, normalize = TRUE)[, 1]

                fc <- bmgarch::forecast(fit_start, ahead = ahead, xC = xC_dat[ i+1, , drop = FALSE], newdata = df_dat[ i+1, ,drop = FALSE])

                loglik[, i+1 ] <-  fc$forecast$log_lik
                approx_elpds_1sap[ i ] <-  .log_sum_exp(lw + loglik[, i+ahead ])
                out[i] <- approx_elpds_1sap[ i ] 
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
            loglik[, i+1 ] <-  fc$forecast$log_lik
        }
        loglik_exact <- loglik[, (L+1):N]
        exact_elpds_1sap <- apply(loglik_exact, 2, .log_mean_exp )
        out <- exact_elpds_1sap
    } else {
        stop("'mode' needs to be either 'forward', 'backward' or 'exact'." )
    }

    refit_info <- cat("Using threshold ", k_thres, 
                      ", model was refit ", length(refits), 
                      " times, at observations", refits)
        
    ## Return relevant objects
    outl <- list( )
    outl$refits <-refits
    outl$ks <- ks
    outl$approx_elpd_1sap <- sum( approx_elpds_1sap, na.rm = TRUE )
    outl$out <- out
    
    outl$loglik <- loglik[, ( L + 1 ):N ]

    return( outl )    
}
