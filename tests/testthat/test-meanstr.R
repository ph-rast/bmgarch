library(bmgarch )
library(testthat )
options(mc.cores=2)


fitVar <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                 parameterization = 'CCC',
                                 meanstructure =  'VAR',
                                 standardize_data = TRUE,
                                 iterations = 10))

test_that("meanstructure returns VAR", {
    expect_equal(fitVar$meanstructure, 2 )
})

fitArma <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                 parameterization = 'CCC',
                                 meanstructure =  'arma',
                                 standardize_data = TRUE,
                                 iterations = 10))

test_that("meanstructure returns VAR", {
    expect_equal(fitArma$meanstructure, 1 )
})

fit <- suppressWarnings( bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
                                 parameterization = 'CCC',
                                 standardize_data = TRUE,
                                 iterations = 10))

test_that("meanstructure returns VAR", {
    expect_equal(fit$meanstructure, 0 )
})

