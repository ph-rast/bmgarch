library(bmgarch )
options(mc.cores=2)

## Fit at least two models to compute model weights
x1 <- x2 <- 0
for( i in 2:100 ) {
    x1[i] <- x1[i-1] + rnorm(1,  0, 1 )
    x2[i] <- x2[i-1] + rnorm(1,  0, 1 )
}

plot(x1)
X <- cbind(x1,  x2 )

fit <- bmgarch(data = X,#stocks[1:300, c("toyota",  "nissan" )],
               Q =  2,
               parameterization = "DCC", standardize_data = FALSE,
               meanstructure = 'VAR',
               iterations = 10000, sampling_algorithm = 'VB')

fit$sampling_algorithm
summary(fit)


plot( rstan::summary(fit$model_fit, pars =  'rts_out')$summary[1:100,'mean'], type =  'l')
lines( x1, col =  'red')
lines(fc$backcast$mean[,,1][,'mean'], col = 'blue')

## same thing:
fitted(fit )[[1]]$mean[,,1]
fc$backcast$mean[,,1][,'mean']


plot(fit, type =  'mean')


fc <- forecast(fit, ahead = 3 )
fc$forecast$mean

plot( fc$backcast$mean[,,1][,'mean'], type =  'l')
lines( x1, col =  'red')

fc$backcast$mean[,,1][,'mean']
x1

fit2 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               P = 2, Q =  2,
               parameterization = "DCC", standardize_data = TRUE,
               iterations = 10000, sampling_algorithm = "VB")
fit2$model_fit[[1]]


## create a bmgarch_list object
blist <- bmgarch_list(fit, fit2 )

## Compute model weights with the default stacking metod
## L is the upper boundary of the time-series before we engage in LFO-CV
mw <- model_weights( blist, L =  98)
print(mw )

fc <- forecast(fit, ahead = 8,
               #xC = cbind(stocks[101:110, c("honda")], stocks[101:110, c("honda")]),
               inc_samples = TRUE,
               newdata = stocks[101:108, c("toyota",  "nissan" )])

plot(fc )

lfob <- loo(fit, mode = 'backward',  L = 99 )
print(lfob)
lfob$out

## obtain model weights for two models
blist <- bmgarch_list(fit, fit2 )
mw <-  model_weights(bmgarch_objects = blist, L =  98)
mw


fc <- forecast( blist, ahead = 1, weights = NULL, xC = cbind(stocks[101, c("honda")], stocks[101, c("honda")]), L =  98)
fc


fcw <- forecast(object = blist, ahead =  5,  weights = mw )
fcw

fc$meta$weights

## Access rstan's model fit object
mf <- fit$model_fit

rstan::check_hmc_diagnostics(mf )
rstan::plot(mf )
