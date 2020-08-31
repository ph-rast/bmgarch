library(bmgarch )
options(mc.cores=2)


## 10 iterations will produce warnings about non-convergence.    
fit <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                 parameterization = "pdBEKK",
                                 standardize_data = TRUE,
                                 iterations = 10))

fit2 <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                  P = 2, Q = 2, 
                                  parameterization = "pdBEKK",
                                  standardize_data = TRUE,
                                  iterations = 10))

source(file =  "../expects.R" , local = TRUE )
