---
title: 'bmgarch: An R-Package for Bayesian Multivariate GARCH models'
tags:
  - R
  - stan
  - GARCH
authors:
  - name: Philippe Rast^[corresponding author]
    orcid: 0000-0003-3630-6629
    affiliation: 1
  - name: Stephen R. Martin
    orcid: 0000-0001-8085-2390
    affiliation: 2
affiliations:
 - name: University of California, Davis
   index: 1
 - name: Comscore, Inc.
   index: 2
date: June 18 2021
bibliography: paper.bib
---

# Introduction
Generalized autoregressive conditional heteroskedasticity models (GARCH) and their multivariate extension (MGARCH) are part of the economists' toolbox ever since their introduction in the early 1980ies [@Bollerslev1986; @Engle1982; @Tsay2013]. Typically, the goal is to generate forecasts of volatility and covolatility for the next day or the near future in time series of assets or market indices. While GARCH models are primarily used in the econometric context, they can be used to capture and forecast heteroskedasticity in any time series. In fact, @Rast2020 presented a parameterization for predicting and capturing within-person variability in human behavior in intensive longitudinal designs.

The main focus of MGARCH models is the estimation of the $p \times p$ conditional covariance matrix $\mathbf{H}_t$ that varies over $t = 1, ... , N$
time points and defines the (co-)volatility of $p$ time series as $\boldsymbol{\epsilon}_t = \mathbf{H}^{1/2}_t \boldsymbol{\eta}_t$, where $\boldsymbol{\epsilon}_t$ is a $p \times 1$ vector of returns or residuals, assuming a process with a conditional mean of zero. $\boldsymbol{\eta}_t$ is a serially independent multivariate white noise process with covariance matrix $\mathbf{I}_p$. While conceptually straightforward, the crux is to define a model that maintains stationarity of the variances and covariances among the time series for all time points as well as positive definiteness for all $\mathbf{H}_t$. As such, a substantial amount of research of the past decades revolved around parameterizations of the conditional covariance matrix that fulfill all those desiderata [for a comparison of some of the most common parameterizations see @DeAlmeida2018]. 

# Statement of need 
While there are a number of readily available packages for univariate GARCH models in R, multivariate implementations are scarcely available. Currently, we are aware of only two packages in the Comprehensive R Archive Network (CRAN). One is `mgarchBEKK` [@mgarchBEKK] which implements BEKK as well as a bivariate asymmetric GARCH model. The other is `rmgarch` [@rmgarch], which includes DCC, GO-GARCH and Copula-GARCH models. Both packages are based on maximum likelihood methods. Moreover, some MGARCH models are implemented in proprietary software (such as Stata), but their focus is mostly limited to CCC, DCC, and VCC models without the option of estimating the location (e.g. time-varying mean returns) and scale (e.g. time-varying volatility) at the same time [@Carnero2014]. Again, these propietary solutions are based on maximum likelihood methods, while in contrast, `bmgarch` implements Bayesian methods for parameter estimation and inference.

At this current time, `bmgarch` implements a CCC [@Bollerslev1990], DCC [@Engle2002], pdBEKK [@Rast2020], and a BEKK [@Engle1995] model based on either Gaussian or Student-t distributed innovations. All parameterizations allow arbitrary ARCH and GARCH orders.
The Bayesian framework is ideally suited to handle the necessary constraints imposed on the MGARCH specifications to both keep the solution stationary and all $t$ conditional covariance matrices $\mathbf{H}_t$ positive definite. Moreover, the model allows one to examine the posterior distribution in order to obtain the full range of Bayesian inference.

The model parameters are estimated using the probabilistic programming language Stan [@stanCppLibrary], which also allows one to use stan related diagnostics implemented in `rstan` [@rstan]. 

# bmgarch

The package is designed to take multivariate time series only. Hence, in terms of data, the minimum requirement is a data frame or matrix with at least two columns, representing the time series. Note that `bmgarch` currently does not support missing values in time series, nor in the predictor variables.

The default behavior in `bmgarch` is to estimate a CCC(1,1) parameterized model of order 1 for the ARCH and GARCH components assuming Student-$t$ distributed innovations and a constant means (zero, typically) structure. The model is fit using 4 chains with 2000 iterations each (half of those are warm-up and the other half are sampling).
As such, the minimum call will be `bmgarch( data = X )` with `X` being a matrix or data frame with at least 2 columns.

In order to include the effect of an external predictor variable on the unconditional and constant covariance term $C$, one can include a predictor `Y` of the same dimension as `X` with `bmgarch( data = X, xC = Y )`. Further, `bmgarch` takes a number of model related arguments governing the type of parameterization (`parameterization = {"CCC", "DCC", "BEKK", "pdBEKK"}`), the order of the GARCH (`P = {1, ..., t-1}`) and ARCH (`Q = {1, ..., t-1}`) processes, the mean structure (`meanstructure = {"constant", "arma"}`) and the type of distribution (`distribution = {"Student_t", "Gaussian"}`). The function also takes arguments with respect to sampling (number of `iterations` and number of `chains`) as well as whether data should be z-standardized before sampling (`standardize_data = {FALSE, TRUE}`) -- this last argument can facilitate the estimation process. 

Objects of the `bmgarch` family can be passed to the `plot()` function, `print()`, and `summary()` for a summary table.

Moreover, each estimated model can be passed to the `forecasting()` function to generate $m$-ahead forecasts of variances, correlations, and, if the ARMA(1,1) is enabled, means. 
Moreover, `bmgarch` integrates model ensemble techniques to generate model weighted forecasts (and estimates) based on Bayesian model averaging and stacking techniques described by @Yao2018 and @Buerkner2020. The forecasted objects can again be plotted and summarized via the corresponding generic functions. 

# Summary
The `bmgarch` package implements four MGARCH parameterizations and estimates all parameters via a Bayesian framework. Sampling from the posterior distribution is based on HMC and is implemented via stan. The `bmgarch` objects can be printed and plotted as well as  passed on to a forecasting function that runs either simple or ensemble forecasts.

# Acknowledgements
The work reported in this publication was supported by the National Institute On Aging of the National Institutes of Health under Award Number R01AG050720 to PR. The content is solely the responsibility 	of the authors and does not necessarily represent the official views of the funding agencies.
