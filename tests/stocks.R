library(bmgarch )
options(mc.cores=2)

## Fit at least two models to compute model weights
x1 <- x2 <- 0
for( i in 2:100 ) {
    x1[i] <- x1[i-1] + rnorm(1,  0, 1 )
    x2[i] <- x2[i-1] + rnorm(1,  0, 1 )
}

devtools::load_all()

.get_target_stan_path( )


out <- standat(data = stocks[1:100, c("toyota",  "nissan" )],
               P = 1, Q = 1, meanstructure = 0, xC = NULL,
               standardize_data = 1, distribution = 0)

stan_data <- out[ c("T", "xC", "rts", "nt", "distribution", "P", "Q", "meanstructure")]
stan_data$rts[100]
stan_data$T

fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan", "honda")],
               Q =  1,
               standardize_data = TRUE,
               parameterization = "BEKK",
               meanstructure = "constant",
               iterations = 50,
               chains = 12,
               sampling_algorithm = 'MCMC')
#               backend =  "cmdstanr")
#               threads =  1,
#               seed =  123,
#               refresh =  0,
#               init =  0,
#               save_latent_dynamics = FALSE,
#               output_dir =  NULL)

## Stuck at this error:
## Error in as.vector(x, "character") : 
##   cannot coerce type 'environment' to vector of type 'character'


summary(fit )
plot(fit )
fc <- forecast(fit,  ahead = 10 )
fc

str(fc$backcast$cor)
fc$backcast$cor[, , "honda_nissan"]

fc$backcast$cor[, "mean", "honda_nissan"]
fc$backcast$cor[, "2.5%", "honda_nissan"]


plot(fc, type = 'var')


fit1 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan", "honda")],
               Q =  1,
               standardize_data = TRUE,
               parameterization = "CCC",
               meanstructure = "constant",
               iterations = 500,
               sampling_algorithm = 'VB',
               backend =  "cmdstanr",
               threads =  1,
               seed =  123,
               refresh =  0,
               init =  0,
               save_latent_dynamics = FALSE,
               output_dir =  NULL)

summary(fit1)

modfits <- bmgarch_list(fit, fit1)

mw <- model_weights(modfits, L = 50, method = 'stacking')
mw

w_fc <- forecast(modfits, ahead = 5, weights = mw )
w_fc
