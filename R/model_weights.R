##' @keywords internal
##' @author philippe
.ll_lfo <- function(x, L = NA, M = M, mode = 'backward') {
    loo( x, type = 'lfo', L = L, M = M, mode = mode)$loglik
}

##' @keywords internal
.rel_eff <- function(ll, x) {
    warmup <- x$model_fit@sim$warmup
    iter <- x$model_fit@sim$iter
    n_chains <- x$model_fit@sim$chains
    loo::relative_eff( exp(ll),
                      chain_id = rep(1:n_chains,  each = iter-warmup ))
}

##' Compute model weights for a list of candidate models based on leave-future-out
##' cross validation (lfocv) expected log-predictive density (elpd).
##' elpd can be approximated via the 'backward' mode described in \insertCite{Buerkner2019;textual}{bmgarch} or via exact cross-validation.    
##' The obtained weights can be passed to the forecast function to obtain weighted forecasts.
##' \code{bmgarch_objects} takes a \code{bmgarch_object} lists. 
##' @title Model weights
##' @param bmgarch_objects list of bmgarch model objects in \code{bmgarch_object}  
##' @param L Minimal length of time series before engaging in lfocv
##' @param M M step head predictions. Defines to what period the LFO-CV should be tuned to. Defaults to M=1.
##' @param method Ensemble methods, 'stacking' (default) or 'pseudobma'
##' @param mode Either 'backward' (default) or 'exact'
##' @return Model weights
##' @details
##' `model_weights()` is a wrapper around the leave-future-out 'lfo' type in `loo.bmgarch()`.
##' The weights can be either obtained from an approximate or exact leave-future-out cross-validation
##' to compute expected log predictive density (ELPD).
##'
##' We can either obtain stacking weights or pseudo-BMA+ weigths as described in \insertCite{Yao2018}{bmgarch}.
##' @examples
##' \dontrun{
##' data(stocks)
##' # Fit at least two models on a subset of the stocks data
##' # to compute model weights
##' fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
##'                parameterization = "DCC", standardize_data = TRUE,
##'                iterations = 500)
##'
##' fit2 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
##'                 P = 2, Q =  2,
##'                parameterization = "DCC", standardize_data = TRUE,
##'                iterations = 500)
##' # create a bmgarch_list object
##' blist <- bmgarch_list(fit, fit2 )
##'
##' # Compute model weights with the default stacking metod
##' # L is the upper boundary of the time-series before we engage in LFO-CV
##' mw <- model_weights( blist, L =  50, method = 'stacking', order = 'backwards' )
##'
##' # Print model weights in the ordert of the bmgarch_list()
##' print(mw)
##' }
##' @references
##'      \insertAllCited{}
##' @export
model_weights <- function(bmgarch_objects = NULL,
                          L = NULL, M = 1,
                          method = "stacking", mode = 'backward') {   

   # if( !is.null(bmgarch_objects ) & !is.null(lfo_objects )) stop( "Supply only 'bmgarch_objects' or 'lfo_objects', not both" )
    if( is.null(bmgarch_objects ) ) stop( "Supply 'bmgarch_objects'" )
    
    if( !is.null(bmgarch_objects ) ) {
    ## Approximates LFO; Ie. results in refitting models.
    ll_list <- lapply(bmgarch_objects, FUN = .ll_lfo, L = L, M = M, mode = mode)
    }
    
    ## Insert lfo_objects
    ## obtain iter, warmup and n_chains from first model
    r_eff_list <- lapply( ll_list, FUN = .rel_eff, bmgarch_objects[[1]] )

    wts <- loo::loo_model_weights( ll_list, method = method,
                                  r_eff_list = r_eff_list,
                                  optim_control = list(reltol=1e-10))
    out <- list()
    out$wts <- wts
    out$ll_list <- ll_list
    out$r_eff_list <- r_eff_list

    attr( out, "class" ) <- "model_weights"
    return(out)
}


##' @title Print method for model_weights
##' @param x Model weights object
##' @param ... Not used.
##' @return model_weights objects with weights, list of log-likelihoods, and r_eff_list 
##' @author philippe
##' @export
print.model_weights <- function(x, ...) {
    print(x$wts)
    return(invisible(x) )
}
