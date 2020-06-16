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

##' .. content for \description{} (no empty lines) ..
##' .. content for \details{} ..
##' @title Model averaging
##' @param bmgarch_objects 
##' @param lfo_objects 
##' @param L 
##' @return Model weights
##' @author philippe
##' @export
model_weights <- function(bmgarch_objects = NULL,  lfo_objects = NULL, L = NULL ) {   

    if( !is.null(bmgarch_objects ) & !is.null(lfo_objects )) stop( "Supply only 'bmgarch_objects' or 'lfo_objects', not both" )
    
    if( !is.null(bmgarch_objects ) ) {
    ## Approximates LFO; Ie. results in refitting models.
    ll_list <- lapply(bmgarch_objects, FUN = .ll_lfo, L = L)
    } else if(!is.null(lfo_objects) ) {
        stop("todo" )
    } else {
        stop("Supply model list for either bmgarch_objects or lfo_objects")
    }
    ## Here, insert lfo_objects
    ## obtain iter, warmup and n_chains from first model
    r_eff_list <- lapply( ll_list, FUN = .rel_eff, bmgarch_objects[[1]] )

    wts <- loo::loo_model_weights( ll_list, method = "stacking",
                                  r_eff_list = r_eff_list,
                                  optim_control = list(reltol=1e-10))
    return(wts)
}
