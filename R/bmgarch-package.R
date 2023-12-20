#' The 'bmgarch' package. 
#' 
#' @description The *bmgarch* package fits Bayesian multivariate GARCH models specified via stan,
#' a C++ package providing HMC methods for full Bayesian inference (cf. [http://mc-stan.org]). The currently implemented parameterizations are DCC(Q,P),
#' CCC(Q,P), and BEKK(Q,P) with arbitrary lags defined in Q, and P. The package provides summaries and plots for the estimates as well
#' as forecasted series with corresponding plots. The fitted objects are rstan class objects that can be inspected and manipulated
#' accordingly.
#'
#' @author Philippe Rast
#' 
#' @docType package
#' @name bmgarch-package
#' @useDynLib bmgarch, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' 
#' 
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#' 
NULL
