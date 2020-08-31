library(bmgarch )
options(mc.cores=2)

## Fit at least two models to compute model weights
fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               parameterization = "DCC", standardize_data = TRUE,
               iterations = 10)

fit2 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               P = 2, Q =  2,
               parameterization = "DCC", standardize_data = TRUE,
               iterations = 10)

## create a bmgarch_list object
blist <- bmgarch_list(fit, fit2 )

## Compute model weights with the default stacking metod
## L is the upper boundary of the time-series before we engage in LFO-CV
mw <- model_weights( blist, L =  98 )
print(mw )

fc <- forecast(fit, ahead = 8,
               #xC = cbind(stocks[101:110, c("honda")], stocks[101:110, c("honda")]),
               inc_samples = TRUE,
               newdata = stocks[101:108, c("toyota",  "nissan" )])

fc

fc$forecast$log_lik

fc <- forecast(fit1, ahead = 8 )

fc

plot(fc )

lfob <- loo(fit, mode = 'backward',  L = 99 )
print(lfob)
lfob$out

## obtain model weights for two models
blist <- bmgarch_list(fit, fit2 )
mw <-  model_weights(bmgarch_objects = blist, L =  98)
mw


fc <- forecast( blist, ahead = 1, weights = NULL, xC = cbind(stocks[101, c("honda")], stocks[101, c("honda")]), L =  80)
fc


fcw <- forecast(object = blist, ahead =  5,  weights = mw )
fcw

fc$meta$weights

## Access rstan's model fit object
mf <- fit$model_fit

rstan::check_hmc_diagnostics(mf )
stan::plot(mf )
