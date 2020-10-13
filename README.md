<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with rmarkdown::render("README.Rmd", output_format = "md_document") -->
<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bmgarch)](https://cran.r-project.org/package=bmgarch)
![build](https://github.com/ph-rast/bmgarch/workflows/build/badge.svg)
![R-CMD-check](https://github.com/ph-rast/bmgarch/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/ph-rast/bmgarch/branch/master/graph/badge.svg)](https://codecov.io/gh/ph-rast/bmgarch)
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
installing the package from CRAN or github. For those who's distro installs
`libnode-dev` instead of `libv8-dev`, run `install.packages("V8")` in R
prior to installing `bmgarch` (during installation`rstan` looks explicitly for V8)

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
    #> Date:  Mon Sep 21 17:55:01 2020
    #> Elapsed time (min):  18.1
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean  sd  mdn  2.5% 97.5%  n_eff Rhat
    #> R_Ng-Ps 0.02 0.5 0.03 -0.92  0.93 833.75 1.01
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.59 0.78 0.28 0.01  2.98 320.54 1.02
    #> var_Ng 1.23 0.40 1.25 0.37  1.98 435.21 1.01
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5%  n_eff Rhat
    #> A_Ps-Ps 0.34 0.10 0.34  0.15  0.55 845.81 1.01
    #> A_Ng-Ps 0.06 0.08 0.07 -0.08  0.20 948.74 1.00
    #> A_Ps-Ng 0.06 0.14 0.06 -0.22  0.32 905.12 1.00
    #> A_Ng-Ng 0.40 0.12 0.40  0.14  0.63 814.42 1.01
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.77 0.21  0.85  0.12  0.94 231.20 1.03
    #> B_Ng-Ps -0.08 0.15 -0.09 -0.41  0.28 259.81 1.00
    #> B_Ps-Ng  0.26 0.37  0.28 -0.68  1.01 444.72 1.00
    #> B_Ng-Ng  0.36 0.20  0.37  0.02  0.73 495.91 1.01
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos -0.02 0.17 -0.01 -0.37  0.30 743.18 1.01
    #> (Intercept)_Neg  0.07 0.12  0.06 -0.14  0.33 444.48 1.01
    #> Phi_Pos-Pos     -0.05 0.36 -0.06 -0.80  0.63 669.45 1.00
    #> Phi_Pos-Neg     -0.11 0.46 -0.13 -0.90  0.74 486.67 1.02
    #> Phi_Neg-Pos     -0.13 0.35 -0.14 -0.80  0.61 370.14 1.01
    #> Phi_Neg-Neg      0.07 0.45  0.09 -0.85  0.85 311.67 1.01
    #> Theta_Pos-Pos   -0.03 0.38 -0.03 -0.74  0.74 649.23 1.00
    #> Theta_Pos-Neg    0.01 0.47  0.04 -0.89  0.82 485.94 1.02
    #> Theta_Neg-Pos    0.16 0.35  0.16 -0.58  0.82 379.07 1.01
    #> Theta_Neg-Neg   -0.08 0.46 -0.09 -0.87  0.84 318.80 1.02
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   52.17   29.16   46.30   13.67  123.43 1375.75    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -802.57    4.75 -802.18 -813.38 -794.72  322.43    1.02

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.41 3.00 -0.34 -6.33  5.21 1916.60    1
    #>    202 -0.16 2.84 -0.21 -5.63  5.44 1752.25    1
    #>    203 -0.10 2.78 -0.09 -5.57  5.37 2019.99    1
    #>    204 -0.13 2.72 -0.10 -5.52  5.20 2037.18    1
    #>    205 -0.12 2.66 -0.18 -5.43  5.27 2032.87    1
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.32 1.50 0.35 -2.61  3.17 1822.75    1
    #>    202 0.27 1.56 0.26 -2.80  3.28 1908.35    1
    #>    203 0.21 1.59 0.20 -2.89  3.30 1965.48    1
    #>    204 0.19 1.63 0.19 -3.15  3.30 1997.93    1
    #>    205 0.17 1.58 0.12 -2.75  3.49 1863.92    1
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 8.15 2.53 7.99 4.06 13.54  582.50    1
    #>    202 7.57 3.51 6.96 3.49 16.40 1093.06    1
    #>    203 7.12 4.51 6.20 3.17 16.47 1544.29    1
    #>    204 6.80 4.97 5.59 2.98 16.92 1934.16    1
    #>    205 6.68 6.67 5.30 2.86 17.82 2025.67    1
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.94 0.39 1.88 1.35  2.85  902.16    1
    #>    202 2.22 0.76 2.04 1.38  4.40 1848.32    1
    #>    203 2.29 1.01 2.08 1.42  4.52 2054.90    1
    #>    204 2.34 1.09 2.07 1.39  5.24 1742.70    1
    #>    205 2.35 1.60 2.01 1.37  5.38 2011.57    1
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.08 0.17 -0.07 -0.43  0.23  620.21    1
    #>    202 -0.09 0.22 -0.09 -0.50  0.37  839.73    1
    #>    203 -0.07 0.22 -0.07 -0.49  0.38 1339.55    1
    #>    204 -0.05 0.22 -0.06 -0.47  0.42 1322.24    1
    #>    205 -0.03 0.21 -0.04 -0.46  0.43 1373.60    1

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
    #> Date:  Mon Sep 21 17:56:08 2020
    #> Elapsed time (min):  0.8
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.11 0.09 0.09 0.00  0.34 2199.63    1
    #> a_h_1,ns   0.08 0.07 0.07 0.00  0.28 2115.85    1
    #> a_h_1,hn   0.11 0.09 0.08 0.00  0.34 2088.59    1
    #> b_h_1,ty   0.46 0.17 0.47 0.10  0.76 1141.02    1
    #> b_h_1,ns   0.37 0.19 0.35 0.06  0.76 1009.63    1
    #> b_h_1,hn   0.39 0.18 0.38 0.08  0.77  994.41    1
    #> c_h_var_ty 0.28 0.12 0.27 0.09  0.56  988.37    1
    #> c_h_var_ns 0.36 0.13 0.36 0.12  0.62 1115.70    1
    #> c_h_var_hn 0.44 0.16 0.43 0.16  0.79  964.89    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.65 0.06 0.65 0.51  0.76 2013.39    1
    #> R_hn-ty 0.73 0.05 0.74 0.62  0.83 1751.63    1
    #> R_hn-ns 0.64 0.07 0.64 0.49  0.76 2171.12    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.08 -0.24  0.07 878.64    1
    #> (Intercept)_nissan -0.01 0.08 -0.01 -0.17  0.15 975.18    1
    #> (Intercept)_honda  -0.02 0.09 -0.02 -0.20  0.16 941.47    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   32.80   23.98   25.76    7.50   97.74 1959.60    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.45    5.31 -177.85 -190.02 -169.51  598.58    1.01

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
    #>    101 0.54 0.11 0.53 0.33  0.77 1364.29    1
    #>    102 0.58 0.16 0.56 0.34  0.90 1826.74    1
    #>    103 0.61 0.21 0.58 0.35  1.05 1937.77    1
    #>    104 0.63 0.22 0.59 0.36  1.22 2067.29    1
    #>    105 0.65 0.33 0.60 0.37  1.22 1958.96    1
    #>    106 0.67 0.44 0.61 0.37  1.31 1820.18    1
    #>    107 0.68 0.33 0.61 0.38  1.38 1699.44    1
    #>    108 0.68 0.39 0.62 0.38  1.41 1930.52    1
    #>    109 0.69 0.44 0.61 0.37  1.51 1668.95    1
    #>    110 0.70 0.49 0.62 0.37  1.53 1746.33    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.61 0.11 0.60 0.42  0.85 2032.13    1
    #>    102 0.65 0.17 0.62 0.42  1.03 1997.43    1
    #>    103 0.66 0.23 0.63 0.42  1.12 1941.16    1
    #>    104 0.68 0.25 0.63 0.42  1.19 2028.20    1
    #>    105 0.68 0.23 0.63 0.42  1.23 1950.30    1
    #>    106 0.68 0.29 0.63 0.43  1.20 1793.60    1
    #>    107 0.68 0.27 0.64 0.43  1.22 1514.76    1
    #>    108 0.68 0.23 0.64 0.43  1.21 1764.11    1
    #>    109 0.68 0.23 0.64 0.42  1.27 1904.80    1
    #>    110 0.69 0.27 0.64 0.43  1.23 1938.49    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.14 0.76 0.53  1.08 1857.11    1
    #>    102 0.84 0.26 0.80 0.53  1.36 1919.94    1
    #>    103 0.87 0.31 0.82 0.53  1.51 2065.57    1
    #>    104 0.90 0.35 0.84 0.54  1.63 2069.19    1
    #>    105 0.92 0.38 0.84 0.54  1.79 2110.90    1
    #>    106 0.92 0.44 0.84 0.55  1.79 2025.72    1
    #>    107 0.94 0.52 0.85 0.55  1.77 1846.19    1
    #>    108 0.94 0.51 0.84 0.55  2.02 2044.99    1
    #>    109 0.93 0.45 0.84 0.55  1.91 2104.66    1
    #>    110 0.93 0.45 0.85 0.55  1.81 1972.44    1

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
    #> Date:  Mon Sep 21 18:10:07 2020
    #> Elapsed time (min):  14.32
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_h_1,ty   0.17 0.12 0.14 0.01  0.46 623.33 1.00
    #> a_h_1,ns   0.10 0.09 0.07 0.00  0.31 395.08 1.01
    #> a_h_1,hn   0.14 0.10 0.12 0.01  0.38 389.50 1.02
    #> b_h_1,ty   0.44 0.16 0.45 0.11  0.73 246.91 1.02
    #> b_h_1,ns   0.42 0.21 0.41 0.07  0.82 204.72 1.01
    #> b_h_1,hn   0.44 0.18 0.43 0.12  0.78 157.04 1.02
    #> c_h_var_ty 0.28 0.12 0.26 0.10  0.54 309.84 1.01
    #> c_h_var_ns 0.31 0.13 0.31 0.08  0.56 230.43 1.01
    #> c_h_var_hn 0.38 0.15 0.37 0.13  0.70 320.45 1.01
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_q 0.21 0.10 0.20 0.04  0.41 402.59    1
    #> b_q 0.23 0.14 0.21 0.01  0.55 479.18    1
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> S_ns-ty 0.61 0.08 0.61 0.43  0.75 358.06 1.01
    #> S_hn-ty 0.73 0.06 0.74 0.60  0.83 512.57 1.01
    #> S_hn-ns 0.62 0.08 0.63 0.46  0.76 325.93 1.00
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota  -0.08 0.08 -0.07 -0.24  0.07 630.52 1.01
    #> (Intercept)_nissan   0.01 0.08  0.02 -0.15  0.17 462.67 1.00
    #> (Intercept)_honda    0.00 0.12  0.00 -0.26  0.23 244.17 1.02
    #> Phi_toyota-toyota    0.03 0.33  0.03 -0.66  0.63 217.46 1.03
    #> Phi_toyota-nissan    0.14 0.43  0.14 -0.70  0.97   9.40 1.17
    #> Phi_toyota-honda     0.09 0.34  0.10 -0.64  0.74 138.40 1.03
    #> Phi_nissan-toyota    0.24 0.38  0.28 -0.53  0.89 165.80 1.03
    #> Phi_nissan-nissan   -0.07 0.41 -0.12 -0.78  0.73  31.31 1.16
    #> Phi_nissan-honda     0.13 0.35  0.15 -0.63  0.74 321.30 1.01
    #> Phi_honda-toyota    -0.19 0.42 -0.22 -0.93  0.60  59.66 1.08
    #> Phi_honda-nissan     0.11 0.42  0.11 -0.70  0.85 187.13 1.07
    #> Phi_honda-honda     -0.18 0.36 -0.20 -0.80  0.56  43.13 1.07
    #> Theta_toyota-toyota -0.12 0.36 -0.13 -0.80  0.63 201.19 1.03
    #> Theta_toyota-nissan  0.01 0.43  0.03 -0.84  0.85  10.68 1.15
    #> Theta_toyota-honda  -0.07 0.34 -0.08 -0.72  0.64 130.49 1.03
    #> Theta_nissan-toyota -0.26 0.38 -0.29 -0.92  0.52 216.91 1.02
    #> Theta_nissan-nissan  0.11 0.38  0.13 -0.64  0.76  43.60 1.12
    #> Theta_nissan-honda  -0.16 0.34 -0.19 -0.80  0.58 346.61 1.00
    #> Theta_honda-toyota  -0.09 0.43 -0.06 -0.91  0.66  38.35 1.09
    #> Theta_honda-nissan   0.02 0.45  0.00 -0.81  0.86 128.56 1.09
    #> Theta_honda-honda    0.31 0.39  0.34 -0.54  0.95  31.63 1.08
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  43.50  27.61  36.83  10.09 112.67 247.91   1.03 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -176.42    5.61 -175.94 -188.58 -166.73  228.88    1.02
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
    #> Using threshold  0.6 , model was refit  3  times, at observations 86 72 59 
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
    #> Using threshold  0.6 , model was refit  3  times, at observations 71 58 50 
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
    #> Using threshold  0.6 , model was refit  13  times, at observations 89 85 81 78 73 69 68 64 61 59 58 57 52

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
    #>    101 -0.09 0.74 -0.08 -1.61  1.36    NA   NA
    #>    102 -0.09 0.74 -0.09 -1.56  1.35    NA   NA
    #>    103 -0.11 0.80 -0.12 -1.64  1.41    NA   NA
    #>    104 -0.09 0.79 -0.09 -1.58  1.45    NA   NA
    #>    105 -0.09 0.82 -0.12 -1.56  1.57    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101  0.02 0.82  0.04 -1.68  1.60    NA   NA
    #>    102 -0.03 0.81 -0.03 -1.61  1.59    NA   NA
    #>    103 -0.02 0.83 -0.03 -1.66  1.66    NA   NA
    #>    104  0.01 0.83  0.02 -1.59  1.64    NA   NA
    #>    105 -0.02 0.81 -0.03 -1.57  1.62    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.02 0.89  0.00 -1.83  1.76    NA   NA
    #>    102 -0.06 0.90 -0.06 -1.87  1.70    NA   NA
    #>    103 -0.03 0.94 -0.02 -1.91  1.90    NA   NA
    #>    104 -0.02 0.93 -0.01 -1.89  1.81    NA   NA
    #>    105 -0.02 0.94 -0.03 -1.95  1.87    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.52 0.10 0.51 0.34  0.74    NA   NA
    #>    102 0.55 0.13 0.53 0.34  0.84    NA   NA
    #>    103 0.60 0.23 0.56 0.35  1.01    NA   NA
    #>    104 0.62 0.50 0.57 0.36  1.10    NA   NA
    #>    105 0.65 0.84 0.59 0.36  1.15    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.63 0.11 0.62 0.46  0.88    NA   NA
    #>    102 0.64 0.13 0.62 0.45  0.97    NA   NA
    #>    103 0.66 0.16 0.64 0.44  1.05    NA   NA
    #>    104 0.67 0.18 0.64 0.44  1.13    NA   NA
    #>    105 0.68 0.19 0.64 0.44  1.15    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.77 0.13 0.76 0.55  1.06    NA   NA
    #>    102 0.79 0.19 0.76 0.52  1.20    NA   NA
    #>    103 0.86 0.29 0.81 0.54  1.52    NA   NA
    #>    104 0.88 0.29 0.82 0.54  1.64    NA   NA
    #>    105 0.90 0.36 0.83 0.55  1.69    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.06 0.65 0.52  0.75    NA   NA
    #>    102 0.64 0.06 0.65 0.52  0.75    NA   NA
    #>    103 0.64 0.06 0.65 0.52  0.75    NA   NA
    #>    104 0.64 0.06 0.65 0.52  0.75    NA   NA
    #>    105 0.64 0.06 0.65 0.52  0.75    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.05 0.74 0.64  0.81    NA   NA
    #>    102 0.73 0.05 0.74 0.64  0.81    NA   NA
    #>    103 0.73 0.05 0.74 0.64  0.81    NA   NA
    #>    104 0.73 0.05 0.74 0.63  0.81    NA   NA
    #>    105 0.73 0.05 0.74 0.63  0.81    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.06 0.64 0.51  0.74    NA   NA
    #>    102 0.64 0.06 0.64 0.51  0.74    NA   NA
    #>    103 0.64 0.06 0.64 0.51  0.74    NA   NA
    #>    104 0.64 0.06 0.64 0.51  0.74    NA   NA
    #>    105 0.64 0.06 0.64 0.51  0.74    NA   NA

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

The preidictors for C are on a log scale in section
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
    #> Date:  Mon Sep 21 20:26:40 2020
    #> Elapsed time (min):  0.58
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.06 0.00  0.35 1872.74 1.00
    #> a_h_1,hn   0.07 0.07 0.05 0.00  0.27 2698.40 1.00
    #> a_h_2,ty   0.09 0.09 0.07 0.00  0.33 2257.48 1.00
    #> a_h_2,hn   0.12 0.12 0.09 0.00  0.46 1563.26 1.00
    #> b_h_1,ty   0.21 0.16 0.17 0.01  0.60 1931.04 1.00
    #> b_h_1,hn   0.18 0.14 0.15 0.01  0.52 1915.13 1.00
    #> b_h_2,ty   0.26 0.17 0.25 0.01  0.61 1477.11 1.01
    #> b_h_2,hn   0.19 0.16 0.15 0.01  0.62 1267.31 1.00
    #> c_h_var_ty 0.22 0.10 0.20 0.07  0.47 1400.22 1.00
    #> c_h_var_hn 0.39 0.16 0.39 0.11  0.73 1242.43 1.00
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.81 2602.17    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.09 -0.09 -0.26  0.09 1558.76    1
    #> (Intercept)_honda  -0.04 0.10 -0.04 -0.23  0.15 1576.46    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.63 0.49 -1.60 -2.68 -0.76 1346.88    1
    #> beta0_hn -1.02 0.46 -0.95 -2.18 -0.31  963.55    1
    #> beta_ty  -0.19 0.40 -0.20 -0.96  0.61 1667.93    1
    #> beta_hn   0.04 0.33  0.07 -0.69  0.66 1315.80    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   45.74   28.46   40.02    8.75  112.63 2897.96    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -129.82    3.99 -129.46 -138.58 -123.09  779.95    1.00

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
    #>    101 0.46 0.11 0.45 0.27  0.73 1555.43    1
    #>    102 0.48 0.16 0.45 0.24  0.82 1771.30    1
    #>    103 0.51 0.24 0.48 0.23  1.01 1905.41    1
    #>    104 0.53 0.28 0.50 0.23  1.02 1872.30    1
    #>    105 0.57 0.36 0.52 0.27  1.23 2005.97    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.15 0.76 0.52  1.08 2144.81    1
    #>    102 0.79 0.24 0.76 0.44  1.36 2319.60    1
    #>    103 0.88 0.35 0.82 0.45  1.72 2211.06    1
    #>    104 0.88 0.35 0.82 0.45  1.74 2101.46    1
    #>    105 0.91 0.39 0.83 0.48  1.86 2088.32    1
