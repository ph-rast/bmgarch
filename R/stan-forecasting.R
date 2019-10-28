##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Forecasting
##' @param object
##' @param ahead
##' @return
##' @author philippe
stan_forecast <- function(object, ahead) {

    ## obtain data from fittd object
    phi0 <- rstan::extract(object$model_fit)[["phi0"]]
    phi <- rstan::extract(object$model_fit)[["phi"]]
    theta <- rstan::extract(object$model_fit)[["theta"]]

    mu <- rstan::extract(object$model_fit)[["mu"]]
    rr <- rstan::extract(object$model_fit)[["rr"]]
    D <- rstan::extract(object$model_fit)[["D"]]
    H <- rstan::extract(object$model_fit)[["H"]]
    
    nu <- rstan::extract(object$model_fit)[["nu"]]

    stan_data_forecast <- list( T = object$TS_length,
                               ahead = ahead,
                               nt = object$nt,
                               Q = object$mgarchQ,
                               P = object$mgarchP,
                               rts = object$RTS_full,
                               xH = object$xH,
                               phi0 = phi0,
                               phi = phi,
                               theta = theta,
                               mu = mu,
                               rr = rr,
                               D = D,
                               H = H,
                               nu = nu,
                               distribution = object$num_dist,
                               meanstructure =  object$meanstructure)

    forecasted <- rstan::sampling(stanmodels$forecast, data = stan_data_forecast,
                                  verbose = TRUE, iter = 100)
    ## Test:
    m <- rstan::stan_model( file = "forecast.stan")
    m
    rstan::sampling(object = m, data =  stan_data_forecast )
}
