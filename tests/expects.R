
test_that("basic model", {
    expect_is( fit, "bmgarch" )                                  
})

test_that("forecast model", {
    fc <- forecast( fit,  ahead = 3 )
    expect_is( fc, "forecast.bmgarch" )                                  
})

test_that("lfo without refits", {
    lfob <- suppressWarnings( loo(fit, mode = 'backward',  L = 99 ) )
    expect_is(lfob, "loo.bmgarch" )
})


test_that("bmgarch model list" , {
    mlist <- bmgarch_list( fit, fit2 )
    expect_is( mlist,"bmgarch_list")
})


test_that("model weights",  {
    mlist <- bmgarch_list( fit, fit2 )
    mw <- suppressWarnings( model_weights( mlist, L = 98 ) )
    expect_is( mw, "model_weights")
})
