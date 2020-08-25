library(bmgarch )
options(mc.cores=2)

fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
#               xC = stocks[1:100, c("honda")],
               meanstructure = "arma",
               parameterization = "CCC", standardize_data = TRUE,
               iterations = 10)

summary(fit)

fit1 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )], 
 #               xC = stocks[1:100, c("honda")],
                meanstructure = "arma",
                parameterization = "DCC", standardize_data = TRUE, #meanstructure = 'arma',
                iterations = 10)

summary(fit1)

fc <- forecast(fit, ahead = 10, xC = cbind(stocks[101:110, c("honda")], stocks[101:110, c("honda")]),
               inc_samples = TRUE, newdata = stocks[101:110, c("toyota",  "nissan" )])

fc
fc$forecast$log_lik

fc <- forecast(fit1, ahead = 10 )

fc

plot(fc )

lfob <- loo(fit, mode = 'backward',  L = 80 )
lfob
lfob$out

## obtain model weights for two models
bmgarch_objects <- bmgarch_list(fit, fit1 )
mw <-  model_weights(bmgarch_objects = bmgarch_objects, L =  90, mode = "forward")
mw
          
fc <- forecast( bmgarch_objects, ahead = 1, weights = NULL, xC = cbind(stocks[101, c("honda")], stocks[101, c("honda")]), L =  80)
fc


fc$meta$weights
