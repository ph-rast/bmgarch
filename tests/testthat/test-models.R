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

        ## test_that("basic model", {
        ##     for(i in 1:4 ) {
        ##      expect_is( fit[[i]], "bmgarch" )   
        ##     }
        ## })

        test_that(paste0("forecast model",i), {
            fc <- forecast( fit,  ahead = 3 )
            expect_is( fc, "forecast.bmgarch" )                                  
        })
        ## test_that("forecast model", {
        ##     for( i in 1:4 ) {
        ##         fc <- forecast( fit[[i]],  ahead = 3 )
        ##         expect_is( fc, "forecast.bmgarch" )
        ##     }
        ## })

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
        
        test_that( paste0(i, "generates non-null output"  ), {
            expect_true( !is.null( summary( fit )$model_summary ) )
        })

        test_that( paste0(i,  "is bmgarch object" ), {
            expect_is( print( fit )$RTS_full, "array" )
        })

        test_that( paste0(i, "generates list of ggplots" ), {
            out <-  plot(fit, askNewPage = FALSE )
            expect_is( out, "list" )    
        })
}

## Simulation of BEKK data:
test_that("Data gets simulated", {
    simdat <- .sim.bekk(N = 100,
                        C = matrix(c(1,0.5,0.5,1), ncol =  2),
                        A = matrix(c(.5,-0.25,0.1,.2), ncol =  2),
                        B = matrix(c(.3,0.1,0.2,.2), ncol =  2),
                        phi = matrix(c(.2,0,0,.2), ncol =  2),
                        theta = matrix(c(.2,0.1,0.1,.2), ncol =  2))
    expect_equal( dim( simdat ) , c(100,  2) )
})


fitvb <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                   sampling_algorithm = 'VB',
                                   parameterization = 'CCC',
                                   standardize_data = TRUE,
                                   iterations = 100))

test_that( "VB returns a VB algo", {
    expect_match( fitvb$sampling_algorithm, "VB"  )
})

fit <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                         parameterization = 'CCC',
                                         standardize_data = TRUE,
                                         iterations = 10))

test_that( "model weights return error if sampling algo does not match", {
    expect_error( model_weights( bmgarch_objects =  bmgarch_list(fit, fitvb), L =  90) )
})
