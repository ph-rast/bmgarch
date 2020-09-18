library(bmgarch )
library(testthat )
options(mc.cores=2)

params <- c('CCC', 'DCC', 'BEKK', 'pdBEKK')

for(i in params ) {
    
    fit <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                     parameterization = i,
                                     standardize_data = TRUE,
                                     iterations = 10))

    fit2 <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                      P =  2, Q =  2,
                                      parameterization = i,
                                      standardize_data = TRUE,
                                      iterations = 10))

    test_that(paste0("basic model ", i), {
        expect_is( fit, "bmgarch" )                                  
    })

    test_that(paste0("forecast model",i), {
        fc <- forecast( fit,  ahead = 3 )
        expect_is( fc, "forecast.bmgarch" )                                  
    })

    test_that(paste0("lfo without refits",i), {
        lfob <- suppressWarnings( loo(fit, mode = 'backward',  L = 99 ) )
        expect_is(lfob, "loo.bmgarch" )
    })

    test_that(paste0("bmgarch model list",i), {
        mlist <- bmgarch_list( fit, fit2 )
        expect_is( mlist,"bmgarch_list")
    })

    test_that(paste0("model weights",i),  {
        mlist <- bmgarch_list( fit, fit2 )
        mw <- suppressWarnings( model_weights( mlist, L = 98 ) )
        expect_is( mw, "model_weights")
    })

}
