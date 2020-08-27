##' @keywords internal
##' @author philippe
.ll_lfo <- function(x, L = NA, mode = 'backward') {
    loo( x, type = 'lfo', L = L, mode = mode)$loglik
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
##' cross validation (lfocv) expected log-predictive predictive density (elpd).
##' elpd can be approximated via the 'backward' mode described in \insertCite{Buerkner2019;textual}{bmgarch} or via exact cross-validation.    
##' The obtained weights can be passed to the forecast function to obtain weighted forecasts.
##' \code{bmgarch_objects} takes a \code{bmgarch_object} lists.
##' @title Model weights
##' @param bmgarch_objects list of bmgarch model objects in \code{bmgarch_object}  
##' @param L Minimal length of time series before engaging in lfocv 
##' @param method Ensemble methods, 'stacking' (default) or 'pseudobma'
##' @param mode Either 'backward' (default) or 'exact'
##' @return Model weights
##' @references
##'      \insertAllCited{}
##' @export
model_weights <- function(bmgarch_objects = NULL,  #lfo_objects = NULL,
                          L = NULL,
                          method = "stacking", mode = 'backward') {   

   # if( !is.null(bmgarch_objects ) & !is.null(lfo_objects )) stop( "Supply only 'bmgarch_objects' or 'lfo_objects', not both" )
    if( is.null(bmgarch_objects ) ) stop( "Supply 'bmgarch_objects'" )
    
    if( !is.null(bmgarch_objects ) ) {
    ## Approximates LFO; Ie. results in refitting models.
    ll_list <- lapply(bmgarch_objects, FUN = .ll_lfo, L = L, mode =  mode)
    }# else if(!is.null(lfo_objects) ) {
        ##if( is.list(lfo_objects ) ) {
        ##    ll_list <- lapply(lfo_objects, FUN = function(x) x$loglik )
        ##}
        ## Need warmup, n_chains and iter form fitted models
        ## if this is to be used, we need to add this to lfo_objects
      # stop("Not yet implemnted" )
   # } else {
    #    stop("Supply model list for either bmgarch_objects or lfo_objects")
    #}

    ## Here, insert lfo_objects
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
