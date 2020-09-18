library(bmgarch )
library(testthat )
options(mc.cores=2)

params <- c('CCC', 'DCC', 'BEKK', 'pdBEKK')

mods <- function(x, P,  Q) {
    bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
            P =  P, Q = Q,
            parameterization = x,
            standardize_data = TRUE,
            iterations = 10)
}

fit <- suppressWarnings( lapply( params, mods, P = 1,  Q = 1 ) )
fit2 <- suppressWarnings( lapply( params, mods, P = 2,  Q = 2 ) )
