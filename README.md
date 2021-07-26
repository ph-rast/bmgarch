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
supports ARMA(1,1) and intercept-only (Constant) mean structures, and a
variety of MGARCH(P,Q) parameterizations. In increasing order of
complexity:

-   CCC(P, Q): Constant Conditional Correlation
-   DCC(P, Q): Dynamic Conditional Correlation
-   BEKK(P, Q): Baba, Engle, Kraft, and Kroner
-   pdBEKK(P, Q): BEKK(P, Q) with positive diagonal constraints

Installation
------------

`bmgarch` is available on CRAN and can be installed with:

    install.packages('bmgarch')

Linux users may need to install `libv8` prior to installing `bmgarch`.
For example, in Ubuntu, run `sudo apt install libv8-dev` before
installing the package from CRAN or github. For those who’s distro
installs `libnode-dev` instead of `libv8-dev`, run
`install.packages("V8")` in R prior to installing `bmgarch` (during
installation`rstan` looks explicitly for V8).

The development version can be installed from
[GitHub](https://github.com/) with:

    devtools::install_github("ph-rast/bmgarch")

How to cite this package
------------------------

Please add following citation when referring to to this package:

Rast, P., Martin, S. R., Liu, S., & Williams, D. R. (in press). A New
Frontier for Studying Within-Person Variability: Bayesian Multivariate
Generalized Autoregressive Conditional Heteroskedasticity Models.
*Psychological Methods*.
<a href="https://doi.org/10.1037/met0000357" class="uri">https://doi.org/10.1037/met0000357</a>

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

    # Fit pdBEKK(1, 1) with ARMA(1,1) on the mean structure.
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
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Mon Jul 12 17:37:47 2021
    #> Elapsed time (min):  14.38
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> R_Ng-Ps 0.01 0.49 0.01 -0.94  0.91 1464.11    1
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.64 0.84 0.29 0.01  3.18 317.46 1.01
    #> var_Ng 1.21 0.42 1.26 0.18  1.91 293.43 1.01
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> A_Ps-Ps 0.34 0.10 0.34  0.13  0.55  797.05    1
    #> A_Ng-Ps 0.06 0.08 0.06 -0.09  0.21 1884.70    1
    #> A_Ps-Ng 0.05 0.14 0.05 -0.25  0.32 2100.56    1
    #> A_Ng-Ng 0.40 0.12 0.40  0.16  0.63  967.36    1
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.75 0.22  0.84  0.12  0.95 300.58 1.01
    #> B_Ng-Ps -0.10 0.15 -0.09 -0.50  0.19 225.38 1.02
    #> B_Ps-Ng  0.29 0.36  0.28 -0.57  1.08 805.29 1.01
    #> B_Ng-Ng  0.36 0.20  0.35  0.03  0.77 705.36 1.00
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_Pos -0.01 0.16 -0.01 -0.37  0.30 1432.45 1.00
    #> (Intercept)_Neg  0.08 0.11  0.07 -0.15  0.31 1712.80 1.00
    #> Phi_Pos-Pos     -0.05 0.36 -0.03 -0.77  0.67  753.00 1.01
    #> Phi_Pos-Neg     -0.12 0.47 -0.12 -0.91  0.76  841.85 1.00
    #> Phi_Neg-Pos     -0.14 0.34 -0.16 -0.76  0.59  803.42 1.00
    #> Phi_Neg-Neg      0.07 0.43  0.08 -0.73  0.82  505.30 1.02
    #> Theta_Pos-Pos   -0.04 0.38 -0.06 -0.78  0.71  732.49 1.00
    #> Theta_Pos-Neg    0.02 0.48  0.04 -0.88  0.83  831.84 1.00
    #> Theta_Neg-Pos    0.16 0.34  0.18 -0.56  0.77  805.93 1.00
    #> Theta_Neg-Neg   -0.08 0.44 -0.09 -0.86  0.74  489.62 1.02
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   52.17   28.15   46.78   13.71  119.87 2897.56    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -802.61    4.60 -802.22 -812.93 -795.00  277.26    1.01

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.49 2.91 -0.54 -6.24  5.17 1975.04    1
    #>    202 -0.25 2.87 -0.20 -6.18  5.39 1903.94    1
    #>    203 -0.07 2.81 -0.08 -5.74  5.39 1868.20    1
    #>    204 -0.08 2.64 -0.09 -5.39  5.08 1848.86    1
    #>    205 -0.05 2.71 -0.07 -5.54  5.33 1817.37    1
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.35 1.48 0.36 -2.44  3.16 1836.81    1
    #>    202 0.24 1.59 0.24 -2.95  3.37 1572.31    1
    #>    203 0.22 1.63 0.21 -2.85  3.39 1983.23    1
    #>    204 0.19 1.63 0.16 -2.94  3.43 1753.18    1
    #>    205 0.10 1.58 0.15 -3.12  3.13 1840.70    1
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 8.06 2.67 7.82 4.07 13.69  637.04    1
    #>    202 7.39 3.49 6.80 3.41 15.57  944.12    1
    #>    203 7.07 4.50 6.04 3.21 16.64 1433.02    1
    #>    204 6.80 5.53 5.51 2.96 18.16 1614.19    1
    #>    205 6.54 5.78 5.19 2.81 17.38 1837.12    1
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.95 0.41 1.89 1.34  2.92 1249.81    1
    #>    202 2.23 0.89 2.02 1.38  4.22 1901.82    1
    #>    203 2.30 1.03 2.07 1.40  4.45 1919.93    1
    #>    204 2.39 2.76 2.05 1.40  5.13 1991.61    1
    #>    205 2.35 1.39 2.04 1.40  5.19 2003.41    1
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.08 0.17 -0.08 -0.43  0.24  818.16    1
    #>    202 -0.08 0.20 -0.08 -0.47  0.35 1199.65    1
    #>    203 -0.07 0.21 -0.07 -0.51  0.38 1395.14    1
    #>    204 -0.05 0.21 -0.06 -0.46  0.41 1487.32    1
    #>    205 -0.05 0.21 -0.05 -0.45  0.38 1616.83    1

    plot(fit.fc, askNewPage = FALSE, type = "var")

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />


    plot(fit.fc, askNewPage = FALSE, type = "cor")

<img src="man/figures/README-forecastPlot-3.png" width="100%" />

Example 2: Stocks
-----------------

Here we use the first 100 days of Stata’s stocks data on daily lagged
returns of three Japanese automakers, Toyota, Nissan, and Honda.

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
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Mon Jul 12 17:38:50 2021
    #> Elapsed time (min):  0.73
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.10 0.09 0.08 0.00  0.34 1878.00    1
    #> a_h_1,ns   0.09 0.08 0.07 0.00  0.26 1262.32    1
    #> a_h_1,hn   0.11 0.09 0.09 0.00  0.31 2288.25    1
    #> b_h_1,ty   0.45 0.18 0.47 0.11  0.76 1117.41    1
    #> b_h_1,ns   0.37 0.19 0.35 0.07  0.78 1097.79    1
    #> b_h_1,hn   0.40 0.18 0.39 0.08  0.75 1283.03    1
    #> c_h_var_ty 0.28 0.12 0.26 0.09  0.57 1036.80    1
    #> c_h_var_ns 0.35 0.13 0.35 0.11  0.61 1105.79    1
    #> c_h_var_hn 0.44 0.16 0.43 0.16  0.77 1264.44    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.64 0.06 0.65 0.51  0.75 1798.57    1
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.82 1785.49    1
    #> R_hn-ns 0.64 0.07 0.64 0.49  0.75 2434.61    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.24  0.06 1167.89    1
    #> (Intercept)_nissan -0.01 0.08 -0.01 -0.16  0.15  996.42    1
    #> (Intercept)_honda  -0.03 0.09 -0.03 -0.20  0.15 1110.12    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   32.37   23.10   25.51    7.10   91.46 2426.88    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.50    5.22 -177.90 -189.77 -169.79  746.39    1.00

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
    #>    101 0.53 0.11 0.52 0.34  0.77 1485.42    1
    #>    102 0.58 0.17 0.56 0.34  0.93 1715.04    1
    #>    103 0.61 0.26 0.58 0.36  1.06 1758.69    1
    #>    104 0.63 0.29 0.59 0.36  1.16 1698.74    1
    #>    105 0.66 0.32 0.60 0.37  1.24 1825.01    1
    #>    106 0.65 0.29 0.60 0.37  1.28 1717.18    1
    #>    107 0.66 0.35 0.60 0.37  1.29 1849.21    1
    #>    108 0.67 0.47 0.60 0.38  1.28 1875.55    1
    #>    109 0.67 0.81 0.61 0.38  1.23 2006.43    1
    #>    110 0.69 0.77 0.60 0.37  1.34 2014.27    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.61 0.11 0.60 0.43  0.85 2178.21    1
    #>    102 0.65 0.16 0.62 0.43  1.04 2025.76    1
    #>    103 0.66 0.21 0.63 0.42  1.14 1963.27    1
    #>    104 0.68 0.24 0.63 0.43  1.16 2130.20    1
    #>    105 0.68 0.23 0.63 0.43  1.16 1994.95    1
    #>    106 0.67 0.25 0.64 0.43  1.12 1225.41    1
    #>    107 0.68 0.30 0.64 0.43  1.16 1668.20    1
    #>    108 0.68 0.28 0.63 0.43  1.18 1795.50    1
    #>    109 0.69 0.50 0.64 0.43  1.20 2002.21    1
    #>    110 0.70 0.84 0.63 0.43  1.26 1703.71    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.14 0.76 0.52  1.07 1598.72    1
    #>    102 0.83 0.22 0.81 0.54  1.31 1916.38    1
    #>    103 0.88 0.32 0.83 0.55  1.61 2037.08    1
    #>    104 0.90 0.32 0.85 0.56  1.62 1971.06    1
    #>    105 0.92 0.45 0.85 0.55  1.81 1659.25    1
    #>    106 0.94 0.74 0.85 0.56  1.78 1492.49    1
    #>    107 0.95 0.82 0.85 0.56  1.96 1696.14    1
    #>    108 0.93 0.46 0.85 0.57  1.71 1587.22    1
    #>    109 0.92 0.36 0.85 0.57  1.71 1747.10    1
    #>    110 0.91 0.34 0.85 0.56  1.76 2091.61    1

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
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Mon Jul 12 17:56:34 2021
    #> Elapsed time (min):  15.91
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_h_1,ty   0.21 0.14 0.19 0.01  0.45   5.40 1.26
    #> a_h_1,ns   0.09 0.07 0.07 0.01  0.28 704.85 1.01
    #> a_h_1,hn   0.16 0.10 0.15 0.01  0.33   5.80 1.25
    #> b_h_1,ty   0.37 0.20 0.38 0.13  0.73   3.81 1.43
    #> b_h_1,ns   0.35 0.21 0.31 0.08  0.80   7.39 1.20
    #> b_h_1,hn   0.38 0.22 0.37 0.11  0.82   4.62 1.33
    #> c_h_var_ty 0.29 0.11 0.31 0.10  0.51  36.01 1.06
    #> c_h_var_ns 0.36 0.14 0.37 0.10  0.55   6.14 1.24
    #> c_h_var_hn 0.40 0.15 0.44 0.12  0.71 144.50 1.03
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> a_q 0.27 0.14 0.26 0.04  0.47  3.10 1.65
    #> b_q 0.20 0.14 0.15 0.01  0.52 11.72 1.12
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> S_ns-ty 0.62 0.08 0.64 0.43  0.74  40.41 1.05
    #> S_hn-ty 0.72 0.06 0.71 0.60  0.84 251.99 1.02
    #> S_hn-ns 0.65 0.08 0.67 0.46  0.76  12.70 1.10
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota  -0.07 0.08 -0.06 -0.24  0.08 307.40 1.03
    #> (Intercept)_nissan   0.04 0.09  0.06 -0.17  0.20  17.39 1.11
    #> (Intercept)_honda   -0.03 0.11 -0.06 -0.23  0.22  48.58 1.07
    #> Phi_toyota-toyota    0.05 0.30  0.17 -0.61  0.59 129.06 1.04
    #> Phi_toyota-nissan   -0.14 0.44 -0.18 -0.69  0.79   5.15 1.33
    #> Phi_toyota-honda     0.00 0.39 -0.03 -0.52  0.77   5.14 1.31
    #> Phi_nissan-toyota    0.34 0.37  0.47 -0.56  0.91  26.05 1.08
    #> Phi_nissan-nissan   -0.16 0.34 -0.23 -0.82  0.58 400.47 1.01
    #> Phi_nissan-honda     0.03 0.38 -0.04 -0.66  0.80  12.54 1.14
    #> Phi_honda-toyota    -0.28 0.36 -0.41 -0.91  0.52 201.73 1.03
    #> Phi_honda-nissan     0.14 0.37  0.26 -0.67  0.84 201.26 1.03
    #> Phi_honda-honda     -0.25 0.41 -0.25 -0.72  0.57   4.43 1.37
    #> Theta_toyota-toyota -0.15 0.33 -0.25 -0.77  0.58 381.46 1.02
    #> Theta_toyota-nissan  0.25 0.40  0.31 -0.62  0.85   8.94 1.20
    #> Theta_toyota-honda   0.07 0.43  0.05 -0.74  0.63   3.51 1.54
    #> Theta_nissan-toyota -0.33 0.36 -0.47 -0.91  0.54 168.22 1.04
    #> Theta_nissan-nissan  0.15 0.32  0.14 -0.54  0.79 481.36 1.00
    #> Theta_nissan-honda  -0.08 0.36 -0.01 -0.83  0.59  16.07 1.11
    #> Theta_honda-toyota   0.05 0.37  0.16 -0.79  0.64  17.25 1.10
    #> Theta_honda-nissan  -0.09 0.42 -0.19 -0.76  0.84  11.80 1.13
    #> Theta_honda-honda    0.40 0.48  0.40 -0.56  0.97   3.96 1.44
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  40.92  26.38  27.91   9.85 109.35  30.73   1.06 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -174.73    6.05 -173.95 -188.23 -167.32    5.56    1.27
    fc <- forecast(fit3, ahead =  10)

    plot( fc,askNewPage = FALSE, type =  'mean' ) 

<img src="man/figures/README-fit3ForecastPlot-1.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-2.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-3.png" width="100%" />

### Compute Model Weights

Obtain model weights with either the stacking or the pseudo BMA method.
These methods are inherited from the `loo` package.

First, gather models to a `bmgarch_list`.

    ## use bmgarch_list function to collect bmgarch objects
    modfits <- bmgarch_list(fit1, fit2, fit3)

Compute model weights with the stacking method (default) and the the
approximate (default) leave-future-out cross validation (LFO CV). `L`
defines the minimal length of the time series before we start engaging
in cross-validation. Eg., for a time series with length 100, `L = 50`
reserves values 51–100 as the cross-validation sample. Note that the
standard is to use the approximate `backward` method to CV as it results
in fewest refits. Exact CV is also available with `exact` but not
encouraged as it results in refitting all CV models.

    mw <- model_weights(modfits, L = 50, method = 'stacking')
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  3  times, at observations 81 72 63 
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  4  times, at observations 81 72 62 58 
    #> 
    #> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  10  times, at observations 89 84 80 74 72 66 60 58 54 52

    ## Return model weights:
    mw
    #> Method: stacking
    #> ------
    #>        weight
    #> model1 0.000 
    #> model2 0.946 
    #> model3 0.054

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
    #>    101 -0.11 0.77 -0.11 -1.65  1.34    NA   NA
    #>    102 -0.11 0.79 -0.11 -1.67  1.44    NA   NA
    #>    103 -0.12 0.79 -0.12 -1.67  1.43    NA   NA
    #>    104 -0.07 0.82 -0.07 -1.71  1.50    NA   NA
    #>    105 -0.08 0.86 -0.07 -1.82  1.58    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.02 0.83 -0.02 -1.67  1.61    NA   NA
    #>    102 -0.01 0.87 -0.03 -1.72  1.69    NA   NA
    #>    103 -0.02 0.86 -0.04 -1.76  1.67    NA   NA
    #>    104  0.02 0.86  0.01 -1.70  1.76    NA   NA
    #>    105  0.04 0.91  0.05 -1.73  1.77    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.05 0.95 -0.05 -2.02  1.86    NA   NA
    #>    102 -0.03 0.95 -0.05 -1.90  1.83    NA   NA
    #>    103 -0.04 0.97 -0.03 -1.95  1.95    NA   NA
    #>    104 -0.06 0.96 -0.04 -1.99  1.81    NA   NA
    #>    105  0.00 1.01  0.01 -2.05  2.04    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.52 0.11 0.52 0.33  0.76    NA   NA
    #>    102 0.55 0.15 0.53 0.33  0.87    NA   NA
    #>    103 0.58 0.17 0.56 0.34  0.99    NA   NA
    #>    104 0.61 0.23 0.57 0.35  1.06    NA   NA
    #>    105 0.62 0.24 0.58 0.35  1.10    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.63 0.11 0.62 0.44  0.88    NA   NA
    #>    102 0.64 0.13 0.62 0.43  0.92    NA   NA
    #>    103 0.66 0.17 0.63 0.43  1.05    NA   NA
    #>    104 0.67 0.19 0.64 0.43  1.08    NA   NA
    #>    105 0.68 0.20 0.64 0.43  1.13    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.78 0.14 0.76 0.54  1.08    NA   NA
    #>    102 0.78 0.18 0.76 0.49  1.19    NA   NA
    #>    103 0.86 0.32 0.80 0.51  1.56    NA   NA
    #>    104 0.87 0.33 0.81 0.51  1.59    NA   NA
    #>    105 0.89 0.38 0.82 0.52  1.77    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.06 0.65 0.51  0.76    NA   NA
    #>    102 0.65 0.06 0.65 0.51  0.76    NA   NA
    #>    103 0.65 0.06 0.65 0.51  0.76    NA   NA
    #>    104 0.65 0.06 0.65 0.51  0.76    NA   NA
    #>    105 0.65 0.06 0.65 0.51  0.76    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.05 0.74 0.62  0.82    NA   NA
    #>    102 0.73 0.05 0.74 0.62  0.82    NA   NA
    #>    103 0.73 0.05 0.74 0.62  0.82    NA   NA
    #>    104 0.73 0.05 0.74 0.62  0.82    NA   NA
    #>    105 0.73 0.05 0.74 0.62  0.82    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.07 0.65 0.49  0.76    NA   NA
    #>    102 0.64 0.07 0.65 0.49  0.76    NA   NA
    #>    103 0.64 0.07 0.65 0.49  0.76    NA   NA
    #>    104 0.64 0.07 0.65 0.49  0.76    NA   NA
    #>    105 0.64 0.07 0.65 0.49  0.76    NA   NA

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

As an example, we will add `nissan` as the predictor for C in a
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

The predictors for C are on a log scale in section
`Exogenous predictor`

    summary(fitx)
    #> Model: CCC-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1
    #> 
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  1000
    #> Chains:  4
    #> Date:  Mon Jul 12 19:55:28 2021
    #> Elapsed time (min):  0.73
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.06 0.00  0.32 1898.77 1.00
    #> a_h_1,hn   0.08 0.08 0.06 0.00  0.30 2245.38 1.00
    #> a_h_2,ty   0.09 0.09 0.07 0.00  0.33 2146.60 1.00
    #> a_h_2,hn   0.12 0.13 0.08 0.00  0.46 1787.97 1.00
    #> b_h_1,ty   0.20 0.16 0.16 0.01  0.59 2384.82 1.00
    #> b_h_1,hn   0.18 0.15 0.14 0.01  0.54 2491.06 1.00
    #> b_h_2,ty   0.26 0.17 0.25 0.01  0.63 1932.13 1.00
    #> b_h_2,hn   0.19 0.16 0.14 0.01  0.63 1398.20 1.00
    #> c_h_var_ty 0.22 0.11 0.21 0.07  0.47 1419.73 1.00
    #> c_h_var_hn 0.39 0.16 0.38 0.12  0.73 1125.11 1.01
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.82 2592.58    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.26  0.07 2108.85    1
    #> (Intercept)_honda  -0.04 0.10 -0.05 -0.23  0.15 1776.44    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.61 0.50 -1.58 -2.69 -0.76 1396.07 1.00
    #> beta0_hn -1.03 0.46 -0.97 -2.12 -0.31  967.57 1.01
    #> beta_ty  -0.19 0.38 -0.19 -0.94  0.59 2510.55 1.00
    #> beta_hn   0.04 0.34  0.06 -0.75  0.67 1860.34 1.00
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   46.11   29.68   39.32    8.81  119.95 3071.53    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -130.12    4.22 -129.76 -139.55 -122.93  631.56    1.00

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
    #>    101 0.47 0.11 0.46 0.28  0.72 1523.69    1
    #>    102 0.48 0.16 0.46 0.24  0.84 2172.20    1
    #>    103 0.51 0.21 0.48 0.23  0.97 2049.67    1
    #>    104 0.54 0.23 0.50 0.25  1.03 1921.62    1
    #>    105 0.58 0.26 0.53 0.27  1.25 2156.70    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.15 0.76 0.52  1.09 2081.04    1
    #>    102 0.80 0.24 0.76 0.43  1.37 2272.83    1
    #>    103 0.88 0.37 0.81 0.45  1.73 2154.06    1
    #>    104 0.89 0.39 0.81 0.45  1.81 1913.36    1
    #>    105 0.92 0.63 0.82 0.48  2.03 2054.87    1

## Community Guidelines

1. Contributions and suggestions to the software are always welcome. Please consult our [contribution guidelines](CONTRIBUTING.md) prior to submitting a pull request.
2. Report issues or problems with the software using github's [issue tracker](https://github.com/ph-rast/bmgarch/issues).
3. Contributors must adhere to the [Code of Conduct](CODE_OF_CONDUCT.md).

## Acknowledgment
This work was supported by the National Institute On Aging of the National Institutes of Health under
Award Number [R01AG050720](https://reporter.nih.gov/project-details/9336761) to PR. The content is solely the responsibility of the authors and does not necessarily represent the official views of the funding agency.
