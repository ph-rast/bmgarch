---
title: 'bmgarch: Bayesian Multivariate GARCH models'
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
    affiliation: 1
affiliations:
 - name: University of California, Davis
   index: 1
date: 8 September 2018
bibliography: paper.bib
---

# Introduction
Generalized autoregressive conditional heteroskedasticity models (GARCH) and their multivariate extension (MGARCH) are part of the economists' tool box ever since their introduction in the early 1980ies [@Bollerslev1986; @Engle1982; @Tsay2013]. Typically, the goal is to generate forecasts of volatility and covolatility for the next day or the near future of assets and markets. While GARCH models are primarily used in the econometric context, they can be used to capture and forecast heteroskedasticity in any time series. In fact, @Rast2020 presented a parameterization  for predicting and capturing within-person variability in human behavior in intensive longitudinal designs. 

The main focus of MGARCH models is the $p \times p$ conditional covariance matrix $\mathbf{H}_t$ that varies over $t = 1, ... , N$
time points and defines the (co-)volatility and of $p$ time series as $\boldsymbol{\epsilon}_t = \mathbf{H}^{1/2}_t \boldsymbol{\eta}_t$. $\boldsymbol{\epsilon}_t$ is a $N \times 1$ vector of returns  or residuals, assuming a process with a conditional mean of zero. $\boldsymbol{\eta}_t$ is a serially independent multivariate white noise process with covariance matrix $\mathbf{I}_p$. While conceptually straightforward, the crux is to define a model that maintains stationarity of the variances and covariances among the time series for all time points as well as positive definiteness for all $\mathbf{H}_t$. As such, a substantial amount of research of the past decades revolves around parameterizations of the conditional covariance matrix that fulfill all those  desiderata [for a comparison of some of the most common parameterizations see @DeAlmeida2018]. 

# Statement of need 
Currently, MGARCH models are mainly implemented in proprietary software, such as Stata, and their focus is mostly limited to CCC, DCC, and VCC models without the option of estimating the location (e.g. time-varying mean returns) and scale (e.g. time-varying volatility) at same time [@Carnero2014]. Moreover, these propietary solutions are based on maximum likelihood methods, while in contrast `bmgarch` implements Bayesian methods to parameter estimation and inference.

At this current time, `bmgarch` implements a CCC [@Bollerslev1990], DCC [@Engle2002], pdBEKK [@Rast2020], and a BEKK [@Engle1995] model based on either Gaussian or Student-t distributed innovations. 
The Bayesian framework is ideally suited to handle the necessary constraints imposed on the MGARCH specifications to both, keep the solutions stationary and all $t$ conditional covariance matrices $\mathbf{H}_t$ positive definite. Moreover, the model allows one to examine the posterior distribution in order obtain the full range of Bayesian inference.

The model parameters are estimated using the probabilistic programming language Stan @stanCppLibrary

# Overview and Getting Started

The default behavior is to assume a constant means term, however, we also offer the estimation of an auto-regressive-moving-average (ARMA 1,1) model of order 1,1. To ensure unbiased estimates when the Student-t distribution is chosen, the ARMA and MGARCH parameters can be estimated simultaneously  @Carnero2014.

Each estimated model can be passed to the forecasting function to generate $m$-ahead forecasts of variances, correlations, and, if ARMA11 is selected, means. 
Moreover, `bmgarch` integrates model ensemble techniques to generate model weighted forecasts (and estimates) based on Bayesian model averaging and stacking techniques described by @Yao2018 and @Buerkner2020.


# Summary

# Acknowledgements
The work reported in this publication was supported by the National Institute On Aging of the National Institutes of Health under Award Number R01AG050720 to PR. The content is solely the responsibility 	of the authors and does not necessarily represent the official views of the funding agencies.
