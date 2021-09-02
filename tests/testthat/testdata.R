library(bmgarch )
library(testthat )
options(mc.cores=2)

context("Data structures of outcomes and predictors" )

test_that("Data vector is not constant", {
    expect_error(
        bmgarch:::standat(data = matrix(c(1, 1, 1, 1, 2, 3 ),  ncol = 2 ),
                         xC = NULL, P = 1,  Q = 1,
                         standardize_data = FALSE, distribution = 1,  meanstructure = 1),
                "Datavector")
})


