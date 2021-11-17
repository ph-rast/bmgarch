<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with rmarkdown::render("README.Rmd", output_format = "md_document") -->
<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bmgarch)](https://cran.r-project.org/package=bmgarch)
![build](https://github.com/ph-rast/bmgarch/workflows/build/badge.svg)
![R-CMD-check](https://github.com/ph-rast/bmgarch/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/ph-rast/bmgarch/branch/master/graph/badge.svg)](https://codecov.io/gh/ph-rast/bmgarch)
[![status](https://joss.theoj.org/papers/fc0db2c4d3bbf1eb18611131f848ed5a/status.svg)](https://joss.theoj.org/papers/fc0db2c4d3bbf1eb18611131f848ed5a)
<!-- badges: end -->

bmgarch
=======

`bmgarch` estimates Bayesian multivariate generalized autoregressive
conditional heteroskedasticity (MGARCH) models. Currently, bmgarch
supports a variety of MGARCH(P,Q) parameterizations and simultaneous
estimation of ARMA(1,1), VAR(1) and intercept-only (Constant) mean
structures. In increasing order of complexity:

-   CCC(P, Q): Constant Conditional Correlation
-   DCC(P, Q): Dynamic Conditional Correlation
-   BEKK(P, Q): Baba, Engle, Kraft, and Kroner
-   pdBEKK(P, Q): BEKK(P, Q) with positive diagonal constraints

Installation
------------

`bmgarch` is available on CRAN and can be installed with:

    install.packages('bmgarch')

### Linux

Linux users may need to install `libv8` prior to installing `bmgarch`.
For example, in Ubuntu, run `sudo apt install libv8-dev` before
installing the package from CRAN or github. For those who’s distro
installs `libnode-dev` instead of `libv8-dev`, run
`install.packages("V8")` in R prior to installing `bmgarch` (during
installation`rstan` looks explicitly for V8).

### Development Version

The development version can be installed from
[GitHub](https://github.com/) with:

    devtools::install_github("ph-rast/bmgarch")

How to cite this package
------------------------

Please add at least one of the following citations when referring to to
this package:

Rast, P., & Martin, S. R. (2021). bmgarch: An R-Package for Bayesian
Multivariate GARCH models. *Journal of Open Source Software*, 6, 3452 -
4354. doi:
<a href="https://doi.org/10.21105/joss.03452" class="uri">https://doi.org/10.21105/joss.03452</a>

Rast, P., Martin, S. R., Liu, S., & Williams, D. R. (in press). A New
Frontier for Studying Within-Person Variability: Bayesian Multivariate
Generalized Autoregressive Conditional Heteroskedasticity Models.
*Psychological Methods*.
<a href="https://doi.org/10.1037/met0000357" class="uri">https://doi.org/10.1037/met0000357</a>;
Preprint-doi:
<a href="https://psyarxiv.com/j57pk/" class="uri">https://psyarxiv.com/j57pk/</a>

Examples:
---------

We present two examples, one with behavioral data and one with stocks
from three major Japanese automakers.

Example 1: Behavioral Data
--------------------------

In this example, we use the pdBEKK(1,1) model for the variances, and an
intercept-only model for the means.

    library(bmgarch)

    data(panas)
    head(panas)
    #>      Pos    Neg
    #> 1 -2.193 -2.419
    #> 2  1.567 -0.360
    #> 3 -0.124 -1.202
    #> 4  0.020 -1.311
    #> 5 -0.150  2.004
    #> 6  3.877  1.008

    ## Fit pdBEKK(1, 1) with ARMA(1,1) on the mean structure.
    fit <- bmgarch(panas,
                   parameterization = "pdBEKK",
                   iterations = 1000,
                   P = 1, Q = 1,
                   distribution = "Student_t",
                   meanstructure = "arma")
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'pdBEKKMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'pdBEKKMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'pdBEKKMGARCH' NOW.

### Parameter estimates

    summary(fit)
    #> Model: pdBEKK-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #> H_t = C + A'[y_(t-1)*y'_(t-1)]A + B'H_(t-1)B
    #> 
    #> Sampling Algorithm:  MCMC
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Wed Nov 17 10:42:49 2021
    #> Elapsed time (min):  19.21
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>          mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #> R_Ng-Ps -0.06 0.38 -0.01 -0.89  0.85 56.84 1.05
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> var_Ps 1.02 0.73 1.26 0.02  2.92 12.69 1.15
    #> var_Ng 1.17 0.33 1.25 0.35  1.83 31.01 1.07
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5% n_eff Rhat
    #> A_Ps-Ps 0.46 0.18 0.42  0.20  0.75  2.27 2.83
    #> A_Ng-Ps 0.05 0.06 0.06 -0.06  0.18  9.52 1.17
    #> A_Ps-Ng 0.10 0.11 0.11 -0.17  0.27  7.81 1.19
    #> A_Ng-Ng 0.39 0.12 0.39  0.17  0.58  4.00 1.42
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.65 0.22  0.71  0.14  0.93   5.27 1.32
    #> B_Ng-Ps -0.05 0.13 -0.03 -0.33  0.24 180.68 1.03
    #> B_Ps-Ng  0.21 0.29  0.27 -0.42  0.93  91.08 1.06
    #> B_Ng-Ng  0.49 0.20  0.61  0.04  0.70   4.05 1.43
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos  0.03 0.13  0.04 -0.33  0.21   7.78 1.19
    #> (Intercept)_Neg  0.07 0.08  0.07 -0.10  0.29 788.15 1.01
    #> Phi_Pos-Pos      0.13 0.34  0.25 -0.77  0.57   5.76 1.27
    #> Phi_Pos-Neg     -0.38 0.42 -0.64 -0.78  0.70   4.69 1.35
    #> Phi_Neg-Pos     -0.21 0.26 -0.16 -0.68  0.43  15.96 1.12
    #> Phi_Neg-Neg      0.26 0.39  0.34 -0.73  0.74   6.16 1.26
    #> Theta_Pos-Pos   -0.29 0.41 -0.42 -0.70  0.75   3.90 1.46
    #> Theta_Pos-Neg    0.34 0.46  0.63 -0.81  0.71   3.79 1.47
    #> Theta_Neg-Pos    0.24 0.28  0.17 -0.43  0.70   9.21 1.18
    #> Theta_Neg-Neg   -0.35 0.44 -0.56 -0.77  0.72   4.49 1.38
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>  mean    sd   mdn  2.5% 97.5% n_eff  Rhat 
    #> 51.82 25.68 45.69 16.84 99.71  4.13  1.41 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -796.43    7.38 -793.28 -811.48 -788.15    2.54    2.21

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -1.16 3.25 -1.15 -7.85  4.89  266.84 1.00
    #>    202 -0.59 3.09 -0.63 -6.50  5.27 1082.14 1.00
    #>    203 -0.44 2.95 -0.50 -6.36  5.49 1533.34 1.00
    #>    204 -0.47 2.77 -0.46 -5.84  4.65 1728.88 1.00
    #>    205 -0.36 2.83 -0.37 -5.69  5.02 2004.18 1.01
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.62 1.51 0.64 -2.20  3.48  121.81 1.01
    #>    202 0.60 1.59 0.61 -2.62  3.66  209.94 1.01
    #>    203 0.52 1.64 0.52 -2.80  3.74  836.84 1.00
    #>    204 0.44 1.70 0.42 -2.88  3.74 1839.11 1.00
    #>    205 0.36 1.67 0.36 -2.86  3.73 1747.95 1.00
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean    sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 9.07  2.59 9.12 4.31 13.45   13.62 1.10
    #>    202 8.21  6.29 6.73 3.52 23.95  442.56 1.03
    #>    203 7.52  8.46 5.81 2.57 22.11 1187.97 1.01
    #>    204 7.17  8.83 5.15 2.33 24.36 1392.74 1.01
    #>    205 7.01 14.12 4.79 2.26 24.62 1908.12 1.01
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.91 0.33 1.89 1.43  2.61   13.10 1.15
    #>    202 2.21 0.76 2.13 1.47  3.96 1598.27 1.00
    #>    203 2.33 0.83 2.16 1.51  4.69 2124.16 1.00
    #>    204 2.44 1.32 2.18 1.48  5.35 1741.05 1.01
    #>    205 2.53 2.73 2.20 1.47  5.68 1991.43 1.00
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #>    201 -0.05 0.14 -0.04 -0.35  0.22 347.38 1.00
    #>    202 -0.02 0.23 -0.04 -0.43  0.56  52.30 1.02
    #>    203  0.01 0.24 -0.02 -0.42  0.63  36.96 1.03
    #>    204  0.03 0.25  0.00 -0.42  0.67  28.53 1.04
    #>    205  0.04 0.25  0.01 -0.40  0.69  24.37 1.05

    plot(fit.fc, askNewPage = FALSE, type = "var")

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />


    plot(fit.fc, askNewPage = FALSE, type = "cor")

<img src="man/figures/README-forecastPlot-3.png" width="100%" />

Example 2: Stocks
-----------------

Here we use the first 100 days (we only base our analyses on 100 days to
reduce wait time – this is not meant to be a serious analysis) of
Stata’s stocks data on daily returns of three Japanese automakers,
Toyota, Nissan, and Honda.

    library(bmgarch)

    data(stocks)
    head(stocks)
    #>         date t       toyota       nissan        honda
    #> 1 2003-01-02 1  0.015167475  0.029470444  0.031610250
    #> 2 2003-01-03 2  0.004820108  0.008173466  0.002679110
    #> 3 2003-01-06 3  0.019958735  0.013064146 -0.001606464
    #> 4 2003-01-07 4 -0.013322592 -0.007444382 -0.011317968
    #> 5 2003-01-08 5 -0.027001143 -0.018856525 -0.016944885
    #> 6 2003-01-09 6  0.011634588  0.016986847  0.013687611

Ease computation by first standardizing the time series

    stocks.z <- scale(stocks[,c("toyota", "nissan", "honda")])
    head(stocks.z )
    #>       toyota     nissan       honda
    #> 1  0.8151655  1.3417896  1.52836901
    #> 2  0.2517820  0.3687089  0.11213515
    #> 3  1.0760354  0.5921691 -0.09765177
    #> 4 -0.7360344 -0.3448866 -0.57304819
    #> 5 -1.4807910 -0.8663191 -0.84849638
    #> 6  0.6228102  0.7714013  0.65102202

    # Fit CCC(1, 1) with constant on the mean structure.
    fit1 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                    parameterization = "CCC",
                    iterations = 1000,
                    P = 1, Q = 1,
                    distribution = "Student_t",
                    meanstructure = "constant")
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.

### Parameter Estimates

    summary( fit1 )
    #> Model: CCC-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1
    #> 
    #> Sampling Algorithm:  MCMC
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Wed Nov 17 10:44:01 2021
    #> Elapsed time (min):  0.9
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.10 0.09 0.08 0.00  0.35 1888.57    1
    #> a_h_1,ns   0.08 0.07 0.06 0.00  0.26 2189.46    1
    #> a_h_1,hn   0.10 0.08 0.09 0.00  0.29 2391.91    1
    #> b_h_1,ty   0.45 0.18 0.46 0.10  0.77 1271.23    1
    #> b_h_1,ns   0.37 0.19 0.35 0.06  0.76 1155.06    1
    #> b_h_1,hn   0.39 0.18 0.38 0.09  0.75 1483.92    1
    #> c_h_var_ty 0.29 0.12 0.27 0.10  0.56 1216.38    1
    #> c_h_var_ns 0.36 0.13 0.36 0.11  0.63 1453.76    1
    #> c_h_var_hn 0.45 0.16 0.43 0.16  0.78 1430.93    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.65 0.06 0.65 0.51  0.75 2407.81    1
    #> R_hn-ty 0.73 0.05 0.74 0.63  0.82 2453.44    1
    #> R_hn-ns 0.64 0.07 0.65 0.50  0.75 2556.38    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.24  0.07 1361.84    1
    #> (Intercept)_nissan -0.01 0.08  0.00 -0.16  0.15 1623.28    1
    #> (Intercept)_honda  -0.02 0.09 -0.02 -0.20  0.17 1510.14    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   32.89   24.58   25.80    7.20   98.90 2315.62    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.38    5.17 -177.93 -189.56 -169.32  726.13    1.00

### Forecasted Values

Forecast volatility 10 days ahead

    fc <- forecast(fit1, ahead = 10 )
    fc
    #> ---
    #> [Variance] Forecast for 10 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.54 0.11 0.53 0.34  0.78 1990.87    1
    #>    102 0.58 0.16 0.56 0.35  0.91 1912.86    1
    #>    103 0.61 0.22 0.58 0.36  1.06 2128.36    1
    #>    104 0.63 0.23 0.59 0.36  1.13 2070.43    1
    #>    105 0.64 0.24 0.60 0.37  1.15 2011.56    1
    #>    106 0.65 0.30 0.60 0.37  1.27 2122.41    1
    #>    107 0.66 0.49 0.61 0.37  1.33 1982.69    1
    #>    108 0.67 0.33 0.61 0.38  1.40 1965.65    1
    #>    109 0.67 0.34 0.61 0.38  1.36 1953.30    1
    #>    110 0.67 0.29 0.62 0.39  1.36 1777.40    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.61 0.11 0.60 0.43  0.84 2199.41    1
    #>    102 0.64 0.16 0.62 0.43  1.01 2165.82    1
    #>    103 0.66 0.19 0.63 0.43  1.15 2258.04    1
    #>    104 0.67 0.19 0.63 0.43  1.16 2218.87    1
    #>    105 0.67 0.20 0.63 0.43  1.16 2160.20    1
    #>    106 0.67 0.21 0.63 0.43  1.20 2097.01    1
    #>    107 0.67 0.23 0.63 0.43  1.17 2093.49    1
    #>    108 0.68 0.27 0.64 0.43  1.20 1662.11    1
    #>    109 0.68 0.27 0.64 0.43  1.26 2134.24    1
    #>    110 0.68 0.24 0.64 0.43  1.26 2095.73    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.14 0.76 0.52  1.07 2113.58    1
    #>    102 0.83 0.22 0.79 0.53  1.36 1926.71    1
    #>    103 0.86 0.27 0.81 0.54  1.45 2023.81    1
    #>    104 0.89 0.34 0.83 0.54  1.61 1781.57    1
    #>    105 0.90 0.35 0.84 0.54  1.67 1921.87    1
    #>    106 0.92 0.43 0.84 0.54  1.74 1773.26    1
    #>    107 0.92 0.47 0.85 0.55  1.68 2081.64    1
    #>    108 0.92 0.39 0.84 0.55  1.80 2181.67    1
    #>    109 0.93 0.46 0.85 0.55  1.78 2044.45    1
    #>    110 0.93 0.39 0.85 0.55  1.78 2023.25    1

    plot(fc,askNewPage = FALSE, type = 'var' )

<img src="man/figures/README-stockForecastPlot-1.png" width="100%" /><img src="man/figures/README-stockForecastPlot-2.png" width="100%" /><img src="man/figures/README-stockForecastPlot-3.png" width="100%" />

### Ensemble Methods

Here we illustrate how to obtain model weights across three models.
These weights will be used to compute weighted forecasts, thus, taking
into account that we do not have a single best model.

Add two additional models, one with CCC(2,2) and a DCC(1,1)

    # Fit CCC(1, 1) with constant on the mean structure.
    fit2 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                    parameterization = "CCC",
                    iterations = 1000,
                    P = 2, Q = 2,
                    distribution = "Student_t",
                    meanstructure = "constant")
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.

    fit3 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                    parameterization = "DCC",
                    iterations = 1000,
                    P = 1, Q = 1,
                    distribution = "Student_t",
                    meanstructure = "arma")
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.

The DCC(1,1) model also incorporates an ARMA(1,1) meanstructure. The
output will have the according information:

    summary( fit3 )
    #> Model: DCC-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]
    #>  R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])
    #> 
    #> Sampling Algorithm:  MCMC
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Wed Nov 17 11:00:01 2021
    #> Elapsed time (min):  14.73
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.17 0.14 0.14 0.01  0.52 1056.07 1.01
    #> a_h_1,ns   0.10 0.09 0.08 0.00  0.34 1189.77 1.00
    #> a_h_1,hn   0.13 0.11 0.11 0.01  0.41 1284.85 1.00
    #> b_h_1,ty   0.44 0.17 0.45 0.11  0.74  906.44 1.00
    #> b_h_1,ns   0.41 0.20 0.39 0.08  0.82  692.83 1.00
    #> b_h_1,hn   0.46 0.19 0.47 0.10  0.83  885.63 1.00
    #> c_h_var_ty 0.28 0.12 0.26 0.10  0.54  935.88 1.00
    #> c_h_var_ns 0.32 0.13 0.32 0.09  0.58  719.80 1.00
    #> c_h_var_hn 0.38 0.16 0.36 0.11  0.72  813.64 1.00
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_q 0.21 0.10 0.20 0.04  0.44 1040.20 1.01
    #> b_q 0.23 0.15 0.21 0.01  0.57  927.87 1.01
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> S_ns-ty 0.60 0.09 0.61 0.40  0.75 1063.10 1.00
    #> S_hn-ty 0.73 0.07 0.74 0.58  0.84  951.12 1.01
    #> S_hn-ns 0.63 0.08 0.63 0.45  0.77 1496.51 1.00
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota  -0.08 0.09 -0.08 -0.26  0.08 792.13 1.00
    #> (Intercept)_nissan   0.01 0.10  0.01 -0.20  0.19 663.45 1.00
    #> (Intercept)_honda   -0.03 0.12 -0.02 -0.27  0.20 644.99 1.00
    #> Phi_toyota-toyota    0.01 0.36  0.02 -0.69  0.68 528.60 1.02
    #> Phi_toyota-nissan    0.02 0.39  0.03 -0.70  0.78 609.47 1.01
    #> Phi_toyota-honda     0.15 0.36  0.16 -0.57  0.87 343.72 1.01
    #> Phi_nissan-toyota    0.27 0.41  0.32 -0.66  0.91 399.28 1.01
    #> Phi_nissan-nissan   -0.15 0.38 -0.17 -0.82  0.64 712.58 1.01
    #> Phi_nissan-honda     0.13 0.41  0.16 -0.76  0.85 446.08 1.01
    #> Phi_honda-toyota    -0.27 0.40 -0.29 -0.94  0.53 558.23 1.01
    #> Phi_honda-nissan     0.14 0.43  0.15 -0.70  0.91 623.26 1.00
    #> Phi_honda-honda     -0.09 0.35 -0.07 -0.76  0.61 665.98 1.00
    #> Theta_toyota-toyota -0.11 0.39 -0.14 -0.83  0.69 416.69 1.02
    #> Theta_toyota-nissan  0.12 0.39  0.13 -0.65  0.82 578.37 1.01
    #> Theta_toyota-honda  -0.13 0.36 -0.14 -0.82  0.57 370.93 1.01
    #> Theta_nissan-toyota -0.27 0.42 -0.34 -0.92  0.71 384.55 1.01
    #> Theta_nissan-nissan  0.16 0.36  0.18 -0.60  0.80 715.11 1.01
    #> Theta_nissan-honda  -0.18 0.40 -0.18 -0.93  0.65 435.38 1.01
    #> Theta_honda-toyota   0.00 0.40  0.00 -0.78  0.73 701.92 1.00
    #> Theta_honda-nissan  -0.02 0.45 -0.03 -0.85  0.86 584.20 1.00
    #> Theta_honda-honda    0.21 0.39  0.21 -0.57  0.91 675.18 1.00
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   44.50   28.32   37.82    9.43  112.97 1642.00    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -177.52    5.99 -177.18 -190.23 -167.03  502.29    1.00
    fc <- forecast(fit3, ahead =  10)

    plot( fc,askNewPage = FALSE, type =  'mean' ) 

<img src="man/figures/README-fit3ForecastPlot-1.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-2.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-3.png" width="100%" />

### Compute Model Weights

Obtain model weights with either the stacking or the pseudo BMA method.
These methods are inherited from the `loo` package.

First, gather models to a `bmgarch_list`.

    ## use bmgarch_list function to collect bmgarch objects
    modfits <- bmgarch_list(fit1, fit2, fit3)

Compute model weights with the stacking method (default) and the
approximate (default) leave-future-out cross validation (LFO CV). `L`
defines the minimal length of the time series before we start engaging
in cross-validation. Eg., for a time series with length 100, `L = 50`
reserves values 51–100 as the cross-validation sample. Note that the
standard is to use the approximate `backward` method to CV as it results
in fewest refits. Exact CV is also available with `exact` but not
encouraged as it results in refitting all CV models.

    mw <- model_weights(modfits, L = 50, method = 'stacking')
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  5  times, at observations 84 77 71 63 51 
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  3  times, at observations 73 65 61 
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'DCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  9  times, at observations 87 84 79 75 74 72 63 60 51

    ## Return model weights:
    mw
    #> Method: stacking
    #> ------
    #>        weight
    #> model1 0.219 
    #> model2 0.781 
    #> model3 0.000

### Weighted Forecasting

Use model weights to obtain weighted forecasts. Here we will forecast 5
days ahead.

    w_fc <- forecast(modfits, ahead = 5, weights = mw )
    w_fc
    #> ---
    #> LFO-weighted forecasts across  3 models.
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.08 0.62 -0.08 -1.29  1.13    NA   NA
    #>    102 -0.10 0.62 -0.09 -1.40  1.10    NA   NA
    #>    103 -0.09 0.66 -0.10 -1.36  1.23    NA   NA
    #>    104 -0.12 0.69 -0.11 -1.56  1.18    NA   NA
    #>    105 -0.10 0.68 -0.11 -1.49  1.20    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101  0.03 0.68  0.04 -1.33  1.38    NA   NA
    #>    102  0.00 0.67  0.01 -1.32  1.34    NA   NA
    #>    103 -0.01 0.70  0.00 -1.34  1.38    NA   NA
    #>    104 -0.03 0.72 -0.02 -1.42  1.39    NA   NA
    #>    105 -0.04 0.72 -0.04 -1.47  1.32    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.02 0.75  0.00 -1.51  1.41    NA   NA
    #>    102 -0.05 0.76 -0.04 -1.59  1.44    NA   NA
    #>    103 -0.02 0.81 -0.01 -1.67  1.51    NA   NA
    #>    104 -0.06 0.83 -0.05 -1.70  1.56    NA   NA
    #>    105 -0.05 0.83 -0.04 -1.71  1.58    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.53 0.09 0.52 0.37  0.72    NA   NA
    #>    102 0.55 0.11 0.54 0.37  0.79    NA   NA
    #>    103 0.59 0.15 0.57 0.39  0.94    NA   NA
    #>    104 0.61 0.17 0.58 0.40  1.00    NA   NA
    #>    105 0.63 0.22 0.59 0.40  1.12    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.63 0.09 0.62 0.47  0.83    NA   NA
    #>    102 0.64 0.11 0.63 0.47  0.87    NA   NA
    #>    103 0.66 0.13 0.64 0.47  0.94    NA   NA
    #>    104 0.66 0.14 0.64 0.46  0.98    NA   NA
    #>    105 0.68 0.16 0.66 0.47  1.01    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.78 0.12 0.77 0.57  1.03    NA   NA
    #>    102 0.79 0.15 0.78 0.55  1.14    NA   NA
    #>    103 0.86 0.22 0.82 0.58  1.43    NA   NA
    #>    104 0.88 0.27 0.83 0.57  1.48    NA   NA
    #>    105 0.91 0.32 0.84 0.59  1.62    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.05 0.65 0.54  0.74    NA   NA
    #>    102 0.65 0.05 0.65 0.54  0.74    NA   NA
    #>    103 0.65 0.05 0.65 0.54  0.74    NA   NA
    #>    104 0.65 0.05 0.65 0.54  0.74    NA   NA
    #>    105 0.65 0.05 0.65 0.54  0.74    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.04 0.74 0.64   0.8    NA   NA
    #>    102 0.73 0.04 0.74 0.64   0.8    NA   NA
    #>    103 0.73 0.04 0.74 0.64   0.8    NA   NA
    #>    104 0.73 0.04 0.74 0.64   0.8    NA   NA
    #>    105 0.73 0.04 0.74 0.64   0.8    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.06 0.64 0.52  0.74    NA   NA
    #>    102 0.64 0.06 0.64 0.52  0.74    NA   NA
    #>    103 0.64 0.06 0.64 0.52  0.74    NA   NA
    #>    104 0.64 0.06 0.64 0.52  0.74    NA   NA
    #>    105 0.64 0.06 0.64 0.52  0.74    NA   NA

Plot the weighted forecast. Save plots into a ggplot object and
post-process

    plt <- plot(w_fc, askNewPage = FALSE, type =  'var' )

<img src="man/figures/README-weightedForecastPlot-1.png" width="100%" /><img src="man/figures/README-weightedForecastPlot-2.png" width="100%" /><img src="man/figures/README-weightedForecastPlot-3.png" width="100%" />


    library( patchwork )
    ( plt$honda  + ggplot2::coord_cartesian(ylim = c(0, 2.5 ) ) ) /
    ( plt$toyota + ggplot2::coord_cartesian(ylim = c(0, 2.5 ) ) ) /
    ( plt$nissan + ggplot2::coord_cartesian(ylim = c(0, 2.5 ) ) ) 
    #> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
    #> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
    #> Coordinate system already present. Adding new coordinate system, which will replace the existing one.

<img src="man/figures/README-weightedForecastPlot-4.png" width="100%" />

### Predictors for Constant Variance (C)

We can add predictors for the constant variance term, c or C, in the
MGARCH model with the option `xC =` The predictors need to be of the
same dimension as the time-series object. For example, with three
time-series of lenght 100, the predictor needs to be entered as a 100 by
3 matrix as well.

To illustrate, we will add `nissan` as the predictor for C in a
bivariate MGARCH:

    # Fit CCC(1, 1) with constant on the mean structure.
    fitx <- bmgarch(stocks.z[1:100, c("toyota", "honda")],
                    xC = stocks.z[1:100, c("nissan", "nissan")],
                    parameterization = "CCC",
                    iterations = 1000,
                    P = 2, Q = 2,
                    distribution = "Student_t",
                    meanstructure = "constant")
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.

The estimates for the predictors for C are on a log scale in section
`Exogenous predictor`:

    summary(fitx)
    #> Model: CCC-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1
    #> 
    #> Sampling Algorithm:  MCMC
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Wed Nov 17 12:45:10 2021
    #> Elapsed time (min):  0.65
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.06 0.00  0.34 1668.39    1
    #> a_h_1,hn   0.08 0.08 0.05 0.00  0.29 1947.72    1
    #> a_h_2,ty   0.10 0.09 0.07 0.00  0.34 1700.01    1
    #> a_h_2,hn   0.12 0.12 0.08 0.00  0.46 1967.79    1
    #> b_h_1,ty   0.20 0.16 0.17 0.01  0.58 2104.55    1
    #> b_h_1,hn   0.18 0.15 0.14 0.01  0.57 1935.18    1
    #> b_h_2,ty   0.26 0.17 0.24 0.01  0.63 1417.28    1
    #> b_h_2,hn   0.19 0.17 0.15 0.01  0.62 1510.62    1
    #> c_h_var_ty 0.22 0.10 0.21 0.07  0.46 1083.03    1
    #> c_h_var_hn 0.40 0.16 0.39 0.12  0.73 1326.13    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.61  0.82 2754.98    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.24  0.07 1594.49    1
    #> (Intercept)_honda  -0.05 0.09 -0.04 -0.23  0.13 1562.95    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.60 0.48 -1.57 -2.60 -0.77 1041.88    1
    #> beta0_hn -1.02 0.47 -0.95 -2.11 -0.32 1062.03    1
    #> beta_ty  -0.20 0.38 -0.20 -0.92  0.58 1729.13    1
    #> beta_hn   0.04 0.32  0.06 -0.66  0.62 1600.12    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   45.76   28.88   39.73    8.86  115.81 2730.22    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -130.18    4.27 -129.61 -139.72 -123.15  674.02    1.01

The predictor results in a linear model (on the log scale) with an
intercept *β*<sub>0</sub> and the effect of the predictor in the slope
*β*<sub>1</sub>.

We can generate forecasts given the known values of the predictor. Note
that the dimension of the predictor needs to match the number of
timepoints that we predict ahead and the number of variables, 5 by 2, in
this example:

    fc2x <- forecast(fitx, ahead = 5, xC = stocks.z[101:105, c("nissan", "nissan")])
    fc2x
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.46 0.11 0.45 0.28  0.69 1719.39    1
    #>    102 0.47 0.16 0.45 0.24  0.82 1923.14    1
    #>    103 0.51 0.23 0.48 0.23  0.96 1826.61    1
    #>    104 0.54 0.29 0.50 0.25  1.07 1918.53    1
    #>    105 0.57 0.27 0.52 0.27  1.18 2027.23    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.15 0.76 0.52  1.11 1591.70    1
    #>    102 0.79 0.24 0.76 0.43  1.35 1851.16    1
    #>    103 0.88 0.34 0.82 0.45  1.82 1929.33    1
    #>    104 0.88 0.36 0.82 0.44  1.72 1992.89    1
    #>    105 0.92 0.51 0.83 0.46  2.10 2073.26    1

### Variational Approximation

The package features the option to use Stan’s variational Bayes
(`sampling_algorithm = "VB"`) algorithm. Currently, this feature is
lagging behind CmdStan’s version and is considered to be experimental
and mostly a placeholder for future improvements.

Community Guidelines
--------------------

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution guidelines](CONTRIBUTING.md) prior
    to submitting a pull request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/ph-rast/bmgarch/issues).
3.  Contributors must adhere to the [Code of
    Conduct](CODE_OF_CONDUCT.md).

Acknowledgment
--------------

This work was supported by the National Institute On Aging of the
National Institutes of Health under Award Number
[R01AG050720](https://reporter.nih.gov/project-details/9336761) to PR.
The content is solely the responsibility of the authors and does not
necessarily represent the official views of the funding agency.
