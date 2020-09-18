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
Generalized autoregressive conditional heteroskedasticity models (GARCH) and their multivariate extension (MGARCH) are part of the economists' tool box ever since their introduction in the early 1980ies [@Bollerslev1986; @Engle1982; @Tsay2013]. The goal is typically to generate forecasts of volatility and covolatility for the next day or the near future of assets and markets. While GARCH models are primarily used in the econometric context, they can be used to capture and forecast heteroskedasticity in any time series. In fact, @Rast2020 presented a parameterization  for predicting and capturing within-person variability in human behavior in intensive longitudinal designs. 

<<<<<<< HEAD
The main focus of MGARCH models is the $p \times p$ conditional covariance matrix $\mathbf{H}_t$ that varies over $t = 1, ... , N$
time points and defines the (co-)volatility and of $p$ time series as $\boldsymbol{\epsilon}_t = \boldsymbol{\epsilon}_t + \mathbf{H}^{1/2}_t \boldsymbol{\eta}_t$. $\boldsymbol{\epsilon}_t$ is a $N \times 1$ vector of returns (for $\boldsymbol{\epsilon} = 0$) or residuals. $\boldsymbol{\eta}_t$ is a serially independent multivariate white noise process with covariance matrix $\mathbf{I}_p$. While conceptually straight forward, the crux is to define a model that maintains stationarity of the variances and covariances among the time series for all time points as well as positive definiteness for all $\mathbf{H}_t$. As such, a substantive amount of research of the past decades revolves around parameterizations of the conditional covariance matrix that fulfill all those  desiderata [for a comparison of some of the most common parameterizations see @DeAlmeida2018]. 
=======
The main focus of MGARCH models is the $p \times p$ conditional covariance matrix $\mathbf{H}_t$ that varies over $t = 1, \ldot , N$ time points and defines the (co-)volatility and of $p$ time series as $\boldsymbol{\epsilon}_t &= \boldsymbol{\epsilon}_t + \mathbf{H}^{1/2}_t \boldsymbol{\eta}_t$. $\boldsymbol{\epsilon}_t$ is a $N \times 1$ vector of returns (for $\boldsymbol{\epsilon} = 0$) or residuals. While conceptually simple, the crux is to define a model that maintains stationarity of the variances and covariances among the time series for all time points as well as positive definiteness for all $\mathbf{H}_t$. As such, a substantive amount of research of the past decades revolves around parameterizations of the conditional covariance matrix that fulfill all those  desiderata [for a comparison of some of the most common parameterizations see @DeAlmeida2018]. 
>>>>>>> 0ae26f456d0b7cbe3ce4eef54dace53b2523c217

# Statement of need 
Currently, MGARCH models are mainly implemented in proprietary software, such as Stata, and their implementation is mostly limited to CCC, DCC, and VCC models without the option of estimating the location (e.g. time-varying mean returns) and scale (e.g. time-varying volatility) at same time [@Carnero2014].
At this time, `bmgarch` implements a CCC [@Bollerslev1990], DCC [@Engle2002], pdBEKK [@Rast2020], and a BEKK [@Engle1995] model with the option of also estimating an auto-regressive-moving-average (ARMA) model of order 1,1 next to a constant means model.
Moreover, `bmgarch` integrates model ensemble techniques to generate model weighted forecasts (and estimates) based on Bayesian model averaging and stacking techniques described by @Yao2018 and @Buerkner2020. 

# Summary
