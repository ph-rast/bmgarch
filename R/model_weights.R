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
model_weights <- function(bmgarch_objects,  lfo_objects = NULL, L = NULL ) {   

   
    ll_list <- lapply(bmgarch_objects, FUN = .ll_lfo, L = L)

    ll_list

    ## obtain iter, warmup and n_chains from first model
    r_eff_list <- lapply( ll_list, FUN = .rel_eff, bmgarch_objects[[1]] )
    r_eff_list
    
#    loo::relative_eff( exp( ll_list[[1]] ), chain_id = rep(1:n_chains,  each = iter-warmup ))

    wts <- loo::loo_model_weights( ll_list, method = "stacking",
                                  r_eff_list = r_eff_list,
                                  optim_control = list(reltol=1e-10))
    return(wts)
}
