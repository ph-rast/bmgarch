library(bmgarch )
options(mc.cores=2)

## Fit at least two models to compute model weights
x1 <- x2 <- 0
for( i in 2:100 ) {
    x1[i] <- x1[i-1] + rnorm(1,  0, 1 )
    x2[i] <- x2[i-1] + rnorm(1,  0, 1 )
}

fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               Q =  1,
               parameterization = "CCC", 
               iterations = 10, sampling_algorithm = 'MCMC')
