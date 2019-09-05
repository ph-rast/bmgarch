#' Standardize input data to facilitate computation
#' 
#' @param data Time-series data
#' @param ahead Number of periods forecasting ahead
#'
#' @return bmgarch object 
#' @export
#' @keywords internal
standat = function(data, ahead, standardize_data, distribution){
    if(dim(data)[1] < dim(data)[2]) data = t(data)
    if(standardize_data) {
    ## Standardize time-series
    stdx = scale(data)
    centered_data = attr(stdx, "scaled:center")
    scaled_data = attr(stdx, "scaled:scale")
    return_standat = list(T = nrow(stdx),
                          rts = stdx,
                          nt = ncol(stdx),
                          ahead = ahead,
                          centered_data = centered_data,
                          scaled_data = scaled_data,
                          distribution = distribution)
    } else {
      ## Unstandardized
      return_standat = list(T = nrow(data),
                            rts = data,
                            nt = ncol(data),
                            ahead = ahead,
                            distribution = distribution)
    }
    return(return_standat)
}

##' Draw samples from a specified multivariate GARCH model, given multivariate time-series.
##'
##' Three paramerization are implemented. The constant conditinal correlation (CCC), the dynamic conditional correlatoin (DCC), and the BEKK.
##' @title Bayesian Multivariate GARCH
##' @param data A time-series or matrix object containing observations at the same interval.
##' @param parameterization A character string specifying the type of of parameterization, must be one of "CCC" (default), "DCC", or "BEKK".
##' @param ahead A number specifying the number of forecasted periods. Defaults to 1.
##' @param iterations A positive integer specifying the number of iterations for each chain (including warmup). The default is 1000
##' @param chains A positive integer specifying the number of Markov chains. The default is 4.
##' @param ... Additional arguments can be ‘chain_id’, ‘init_r’, ‘test_grad’, ‘append_samples’, ‘refresh’, ‘enable_random_init’. See the documentation in ‘stan’.
##' @return An object of S4 class ‘stanfit’ representing the fitted results.
##' @author philippe
##' @export
bmgarch = function(data, parameterization = 'CCC', ahead = 1, iterations = 1000, chains = 4, standardize_data = TRUE,
                   distribution = 'Normal', ...) {
    num_dist = NA
    if ( distribution == 'Normal' ) num_dist = 0 else {
               if ( distribution == 'student_t' ) num_dist = 1 else warning( '\n\n Specify distribution: Normal or student_t \n\n', immediate. = TRUE) }
    return_standat = standat(data, ahead, standardize_data, distribution = num_dist )
    stan_data  = return_standat[ c('T', 'rts', 'nt', 'ahead', 'distribution')]

  if(parameterization == 'CCC') model_fit <- rstan::sampling(stanmodels$CCCMGARCH, data = stan_data,
                                                      verbose = TRUE,
                                                      iter = iterations,
                                                      control = list(adapt_delta = .99),
                                                      init_r = 1,
                                                      chains = chains) else {
  if( parameterization == 'DCC' ) model_fit <- rstan::sampling(stanmodels$DCCMGARCH, data = stan_data,
                                                      verbose = TRUE,
                                                      iter = iterations,
                                                      control = list(adapt_delta = .99),
                                                      init_r = 1,
                                                      chains = chains) else {
  if( parameterization == 'BEKK' ) model_fit <- rstan::sampling(stanmodels$BEKKMGARCH, 
                                                      data = stan_data,
                                                      verbose = TRUE,
                                                      iter = iterations,
                                                      control = list(adapt_delta = .99),
                                                      init_r = 1,
                                                      chains = chains) else {
  warning( 'Not a valid model specification. Select CCC, DCC, or BEKK.' )}
                                                                       }
                                                                       }
    ## Model fit is based on standardized values.
    mns = return_standat$centered_data
    sds = return_standat$scaled_data
    ## Values could be converted to original scale using something like this on the estimates
    ## orig_sd = stan_data$rts %*% diag(sds)
    ## orig_scale = orig_sd + array(rep(mns, each = aussi[[1]]$T), dim = c(aussi[[1]]$T, aussi[[1]]$nt) )
    return_fit <- list(model_fit = model_fit,
                       param = parameterization,
                       distribution = distribution,
                       iter = iterations,
                       chains = chains,
                       elapsed_time = rstan::get_elapsed_time(model_fit),
                       date = date(),
                       ahead = ahead,
                       nt = stan_data$nt,
                       TS_length = stan_data$T,
                       TS_names = colnames(stan_data$rts),
                       RTS_last = stan_data$rts[stan_data$T,],
                       RTS_full = stan_data$rts)
    class(return_fit) <- "bmgarch"
    return(return_fit)
}
