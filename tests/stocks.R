library(bmgarch )
options(mc.cores=2)

fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               parameterization = "DCC", standardize_data = TRUE,
               iterations = 10)

summary(fit)

fit1 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               P = 2, Q =  2,
               parameterization = "DCC", standardize_data = TRUE,
               iterations = 10)

summary(fit1)

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
lfob
lfob$out

## obtain model weights for two models
blist <- bmgarch_list(fit, fit1 )
blist
mw <-  model_weights(bmgarch_objects = blist, L =  98)
mw


fc <- forecast( bmgarch_objects, ahead = 1, weights = NULL, xC = cbind(stocks[101, c("honda")], stocks[101, c("honda")]), L =  80)
fc


fc$meta$weights
