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

Note that during installation all models are compiled – this might take
a while.

The development version can be installed from
[GitHub](https://github.com/) with:

    devtools::install_github("ph-rast/bmgarch")

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
    #> Date:  Fri Sep 18 13:05:05 2020
    #> Elapsed time (min):  15.99
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean   sd  mdn  2.5% 97.5%  n_eff Rhat
    #> R_Ng-Ps 0.01 0.49 0.01 -0.91  0.91 818.31 1.01
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.58 0.79 0.26 0.01  3.08 281.46 1.01
    #> var_Ng 1.24 0.39 1.27 0.39  1.96 514.50 1.00
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> A_Ps-Ps 0.34 0.10 0.34  0.14  0.55 1187.14 1.00
    #> A_Ng-Ps 0.06 0.07 0.06 -0.08  0.20  910.18 1.00
    #> A_Ps-Ng 0.06 0.14 0.06 -0.23  0.33  752.67 1.01
    #> A_Ng-Ng 0.40 0.12 0.40  0.15  0.62  865.49 1.01
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.77 0.21  0.85  0.10  0.95 169.72 1.03
    #> B_Ng-Ps -0.08 0.14 -0.09 -0.39  0.27 273.32 1.02
    #> B_Ps-Ng  0.25 0.37  0.27 -0.78  0.98 293.00 1.01
    #> B_Ng-Ng  0.36 0.20  0.35  0.03  0.74 718.08 1.00
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos -0.01 0.16 -0.01 -0.35  0.31 632.59 1.01
    #> (Intercept)_Neg  0.08 0.12  0.07 -0.16  0.33 727.60 1.00
    #> Phi_Pos-Pos     -0.04 0.37 -0.03 -0.83  0.65 386.74 1.01
    #> Phi_Pos-Neg     -0.11 0.47 -0.12 -0.92  0.78 377.71 1.01
    #> Phi_Neg-Pos     -0.16 0.34 -0.18 -0.77  0.55 507.04 1.01
    #> Phi_Neg-Neg      0.04 0.43  0.04 -0.78  0.80 436.46 1.01
    #> Theta_Pos-Pos   -0.05 0.39 -0.06 -0.77  0.77 382.33 1.01
    #> Theta_Pos-Neg    0.02 0.48  0.04 -0.90  0.83 371.64 1.01
    #> Theta_Neg-Pos    0.18 0.34  0.20 -0.53  0.80 476.18 1.01
    #> Theta_Neg-Neg   -0.05 0.44 -0.04 -0.85  0.79 417.52 1.01
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   53.63   27.98   48.82   14.05  120.10 1356.86    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -802.62    5.10 -801.99 -814.65 -794.18  212.97    1.04

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.43 3.03 -0.49 -6.27  5.49 1769.15    1
    #>    202 -0.18 2.91 -0.25 -6.01  5.92 2147.19    1
    #>    203  0.03 2.73  0.05 -5.57  5.09 1877.05    1
    #>    204 -0.17 2.63 -0.22 -5.15  5.45 1902.87    1
    #>    205  0.01 2.69  0.04 -5.04  5.35 2008.84    1
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.36 1.45 0.32 -2.48  3.25 1726.54    1
    #>    202 0.27 1.60 0.28 -2.93  3.34 1777.92    1
    #>    203 0.18 1.62 0.15 -3.00  3.52 2089.03    1
    #>    204 0.22 1.57 0.23 -2.97  3.31 1530.69    1
    #>    205 0.07 1.55 0.10 -2.96  3.13 1858.73    1
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 8.22 2.60 8.10 4.03 13.79  554.01    1
    #>    202 7.67 3.48 7.10 3.43 16.29  961.67    1
    #>    203 7.28 4.09 6.32 3.25 17.64 1226.97    1
    #>    204 6.93 4.39 5.80 3.08 18.10 1379.45    1
    #>    205 6.60 4.41 5.39 2.93 18.57 1585.14    1
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.94 0.40 1.89 1.34  2.92  965.77    1
    #>    202 2.23 0.81 2.04 1.38  4.25 1555.31    1
    #>    203 2.32 1.06 2.07 1.41  4.74 1750.87    1
    #>    204 2.37 1.14 2.08 1.40  5.27 1910.88    1
    #>    205 2.34 1.27 2.04 1.39  5.35 1618.09    1
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.08 0.17 -0.08 -0.42  0.25  926.54    1
    #>    202 -0.08 0.22 -0.09 -0.49  0.39 1100.38    1
    #>    203 -0.07 0.22 -0.07 -0.48  0.39 1327.51    1
    #>    204 -0.05 0.22 -0.05 -0.47  0.40 1395.87    1
    #>    205 -0.04 0.21 -0.05 -0.45  0.40 1639.85    1

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
    #> Date:  Fri Sep 18 13:06:23 2020
    #> Elapsed time (min):  0.96
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.11 0.10 0.08 0.00  0.37 1956.97    1
    #> a_h_1,ns   0.09 0.08 0.07 0.00  0.28 3208.74    1
    #> a_h_1,hn   0.11 0.08 0.10 0.01  0.32 3292.94    1
    #> b_h_1,ty   0.45 0.17 0.46 0.10  0.76 1942.82    1
    #> b_h_1,ns   0.37 0.19 0.36 0.06  0.76 1423.57    1
    #> b_h_1,hn   0.38 0.17 0.38 0.07  0.73 1838.74    1
    #> c_h_var_ty 0.29 0.12 0.27 0.09  0.56 1692.58    1
    #> c_h_var_ns 0.36 0.13 0.36 0.12  0.63 1388.86    1
    #> c_h_var_hn 0.45 0.16 0.44 0.17  0.77 1934.19    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.65 0.07 0.65 0.51  0.76 2417.99    1
    #> R_hn-ty 0.73 0.05 0.74 0.63  0.82 2381.08    1
    #> R_hn-ns 0.64 0.07 0.65 0.50  0.76 2570.81    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.25  0.07 1350.33    1
    #> (Intercept)_nissan -0.01 0.09 -0.01 -0.18  0.16 1432.93    1
    #> (Intercept)_honda  -0.02 0.10 -0.03 -0.21  0.17 1220.80    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   32.92   24.49   24.20    7.28   95.34 2114.05    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.38    5.08 -178.01 -189.62 -169.63  703.05    1.01

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
    #>    101 0.54 0.11 0.53 0.33  0.78 1846.54    1
    #>    102 0.59 0.19 0.57 0.35  1.00 1939.94    1
    #>    103 0.62 0.22 0.59 0.35  1.09 1827.29    1
    #>    104 0.63 0.25 0.60 0.36  1.16 1924.47    1
    #>    105 0.64 0.23 0.60 0.37  1.22 1977.78    1
    #>    106 0.66 0.41 0.60 0.38  1.31 2003.62    1
    #>    107 0.68 0.59 0.61 0.37  1.38 1931.29    1
    #>    108 0.68 0.44 0.60 0.38  1.51 1851.03    1
    #>    109 0.68 0.51 0.61 0.38  1.44 1819.30    1
    #>    110 0.69 0.52 0.61 0.38  1.38 1847.50    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.62 0.11 0.61 0.43  0.86 1785.54    1
    #>    102 0.65 0.20 0.63 0.42  1.04 1782.78    1
    #>    103 0.67 0.25 0.64 0.42  1.14 1875.90    1
    #>    104 0.69 0.30 0.64 0.42  1.19 1978.31    1
    #>    105 0.69 0.26 0.64 0.43  1.22 2004.39    1
    #>    106 0.69 0.25 0.64 0.43  1.29 1905.76    1
    #>    107 0.68 0.24 0.64 0.43  1.24 1986.51    1
    #>    108 0.69 0.25 0.64 0.43  1.24 1930.03    1
    #>    109 0.68 0.23 0.64 0.43  1.20 1809.34    1
    #>    110 0.69 0.28 0.64 0.43  1.24 1951.77    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.76 0.14 0.75 0.53  1.06 2494.87    1
    #>    102 0.85 0.26 0.80 0.54  1.51 2032.03    1
    #>    103 0.88 0.31 0.82 0.55  1.58 2068.14    1
    #>    104 0.90 0.35 0.84 0.56  1.77 2147.53    1
    #>    105 0.92 0.40 0.84 0.56  1.81 1839.33    1
    #>    106 0.93 0.55 0.85 0.56  1.76 2106.19    1
    #>    107 0.94 0.63 0.85 0.55  1.89 1878.88    1
    #>    108 0.93 0.43 0.85 0.56  1.92 2017.77    1
    #>    109 0.93 0.43 0.85 0.57  1.86 2107.94    1
    #>    110 0.94 0.41 0.86 0.56  1.91 2033.04    1

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
    #> Date:  Fri Sep 18 13:22:03 2020
    #> Elapsed time (min):  14.57
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.17 0.14 0.14 0.01  0.55 1209.50    1
    #> a_h_1,ns   0.10 0.08 0.08 0.00  0.30 1980.04    1
    #> a_h_1,hn   0.13 0.10 0.11 0.01  0.37 1359.36    1
    #> b_h_1,ty   0.45 0.17 0.46 0.11  0.75 1385.97    1
    #> b_h_1,ns   0.42 0.21 0.40 0.08  0.83 1242.67    1
    #> b_h_1,hn   0.45 0.19 0.46 0.11  0.80 1375.71    1
    #> c_h_var_ty 0.27 0.12 0.25 0.09  0.53 1330.93    1
    #> c_h_var_ns 0.32 0.13 0.31 0.08  0.58 1426.16    1
    #> c_h_var_hn 0.39 0.16 0.38 0.13  0.73 1376.59    1
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_q 0.21 0.10 0.19 0.04  0.43 1334.58    1
    #> b_q 0.24 0.15 0.22 0.02  0.56 1532.29    1
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> S_ns-ty 0.60 0.09 0.61 0.41  0.75 1758.11    1
    #> S_hn-ty 0.73 0.07 0.74 0.58  0.84 2003.39    1
    #> S_hn-ns 0.63 0.08 0.64 0.44  0.77 1383.31    1
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota  -0.07 0.09 -0.07 -0.26  0.10 1022.94 1.00
    #> (Intercept)_nissan   0.02 0.09  0.01 -0.16  0.21  997.81 1.00
    #> (Intercept)_honda   -0.01 0.11 -0.01 -0.24  0.22 1197.43 1.00
    #> Phi_toyota-toyota    0.03 0.33  0.03 -0.62  0.66  707.25 1.00
    #> Phi_toyota-nissan    0.01 0.39  0.00 -0.73  0.77  534.43 1.01
    #> Phi_toyota-honda     0.10 0.36  0.12 -0.67  0.76  765.92 1.00
    #> Phi_nissan-toyota    0.25 0.40  0.28 -0.59  0.92  546.62 1.01
    #> Phi_nissan-nissan   -0.14 0.40 -0.16 -0.85  0.65  733.75 1.00
    #> Phi_nissan-honda     0.10 0.42  0.13 -0.76  0.84  472.03 1.01
    #> Phi_honda-toyota    -0.25 0.40 -0.27 -0.94  0.58  777.45 1.00
    #> Phi_honda-nissan     0.15 0.42  0.14 -0.66  0.91  704.52 1.01
    #> Phi_honda-honda     -0.11 0.33 -0.11 -0.74  0.54  923.12 1.01
    #> Theta_toyota-toyota -0.12 0.36 -0.13 -0.76  0.59  623.31 1.01
    #> Theta_toyota-nissan  0.14 0.39  0.16 -0.66  0.86  513.71 1.01
    #> Theta_toyota-honda  -0.09 0.35 -0.10 -0.75  0.63  796.41 1.00
    #> Theta_nissan-toyota -0.26 0.41 -0.30 -0.91  0.62  530.49 1.01
    #> Theta_nissan-nissan  0.16 0.37  0.17 -0.57  0.83  729.32 1.00
    #> Theta_nissan-honda  -0.14 0.40 -0.15 -0.90  0.68  542.98 1.01
    #> Theta_honda-toyota  -0.02 0.40  0.00 -0.80  0.71  786.10 1.01
    #> Theta_honda-nissan  -0.04 0.45 -0.03 -0.86  0.84  622.78 1.01
    #> Theta_honda-honda    0.22 0.37  0.23 -0.51  0.91  817.25 1.01
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   44.88   29.70   37.74    8.09  117.22 2148.30    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -177.22    5.64 -176.92 -188.65 -166.92  656.22    1.00
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
    #> Using threshold  0.6 , model was refit  3  times, at observations 87 72 61 
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
    #> Using threshold  0.6 , model was refit  3  times, at observations 81 72 61 
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
    #> Using threshold  0.6 , model was refit  12  times, at observations 89 85 81 74 73 68 66 64 61 58 57 53

    ## Return model weights:
    mw
    #> Method: stacking
    #> ------
    #>        weight
    #> model1 0.000 
    #> model2 0.897 
    #> model3 0.103

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
    #>    101 -0.11 0.69 -0.12 -1.46  1.31    NA   NA
    #>    102 -0.09 0.70 -0.11 -1.49  1.32    NA   NA
    #>    103 -0.08 0.74 -0.08 -1.55  1.36    NA   NA
    #>    104 -0.10 0.74 -0.11 -1.53  1.36    NA   NA
    #>    105 -0.08 0.79 -0.07 -1.70  1.46    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.01 0.78 -0.01 -1.55  1.51    NA   NA
    #>    102  0.00 0.76 -0.01 -1.55  1.55    NA   NA
    #>    103 -0.01 0.79  0.00 -1.57  1.52    NA   NA
    #>    104 -0.01 0.78  0.01 -1.58  1.49    NA   NA
    #>    105  0.01 0.80 -0.01 -1.51  1.62    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.05 0.84 -0.07 -1.71  1.60    NA   NA
    #>    102 -0.02 0.84 -0.02 -1.71  1.65    NA   NA
    #>    103 -0.03 0.88 -0.02 -1.79  1.69    NA   NA
    #>    104 -0.02 0.89 -0.01 -1.77  1.75    NA   NA
    #>    105 -0.01 0.93  0.00 -1.87  1.84    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.52 0.10 0.51 0.34  0.74    NA   NA
    #>    102 0.55 0.13 0.54 0.35  0.84    NA   NA
    #>    103 0.59 0.17 0.56 0.36  1.02    NA   NA
    #>    104 0.61 0.20 0.58 0.37  1.06    NA   NA
    #>    105 0.63 0.23 0.59 0.37  1.17    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.63 0.11 0.62 0.45  0.87    NA   NA
    #>    102 0.64 0.14 0.62 0.45  0.95    NA   NA
    #>    103 0.66 0.17 0.64 0.44  1.04    NA   NA
    #>    104 0.67 0.19 0.64 0.44  1.08    NA   NA
    #>    105 0.68 0.20 0.64 0.45  1.13    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.77 0.13 0.76 0.56  1.07    NA   NA
    #>    102 0.78 0.17 0.76 0.52  1.13    NA   NA
    #>    103 0.85 0.25 0.80 0.54  1.41    NA   NA
    #>    104 0.86 0.26 0.81 0.54  1.51    NA   NA
    #>    105 0.89 0.30 0.83 0.55  1.61    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.06 0.65 0.52  0.75    NA   NA
    #>    102 0.64 0.06 0.65 0.51  0.74    NA   NA
    #>    103 0.64 0.06 0.64 0.51  0.74    NA   NA
    #>    104 0.64 0.06 0.64 0.51  0.75    NA   NA
    #>    105 0.64 0.06 0.65 0.51  0.74    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.05 0.73 0.62  0.81    NA   NA
    #>    102 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    103 0.73 0.05 0.73 0.63  0.82    NA   NA
    #>    104 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    105 0.73 0.05 0.73 0.62  0.82    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.06 0.64 0.51  0.75    NA   NA
    #>    102 0.64 0.06 0.64 0.50  0.75    NA   NA
    #>    103 0.64 0.06 0.64 0.50  0.75    NA   NA
    #>    104 0.64 0.06 0.64 0.50  0.75    NA   NA
    #>    105 0.64 0.06 0.64 0.50  0.75    NA   NA

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
    #> Date:  Fri Sep 18 15:32:27 2020
    #> Elapsed time (min):  0.66
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.06 0.00  0.34 1701.40    1
    #> a_h_1,hn   0.08 0.08 0.05 0.00  0.29 2281.21    1
    #> a_h_2,ty   0.10 0.09 0.07 0.00  0.33 2011.41    1
    #> a_h_2,hn   0.12 0.12 0.08 0.00  0.47 1377.75    1
    #> b_h_1,ty   0.21 0.17 0.17 0.01  0.60 1807.46    1
    #> b_h_1,hn   0.18 0.15 0.15 0.01  0.55 2065.11    1
    #> b_h_2,ty   0.26 0.17 0.23 0.01  0.63 2109.38    1
    #> b_h_2,hn   0.19 0.17 0.14 0.00  0.61 1150.76    1
    #> c_h_var_ty 0.22 0.10 0.21 0.07  0.46 1326.22    1
    #> c_h_var_hn 0.39 0.16 0.38 0.11  0.74  829.50    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.81 2922.53    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.25  0.08 1568.29    1
    #> (Intercept)_honda  -0.05 0.09 -0.05 -0.22  0.14 1557.55    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.61 0.48 -1.58 -2.67 -0.78 1291.03    1
    #> beta0_hn -1.04 0.48 -0.97 -2.21 -0.29  651.36    1
    #> beta_ty  -0.19 0.38 -0.19 -0.96  0.56 1851.72    1
    #> beta_hn   0.03 0.32  0.06 -0.68  0.60 1221.99    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   45.69   30.06   38.48    7.80  120.25 2721.89    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -129.85    4.01 -129.55 -138.39 -122.82  785.98    1.00

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
    #>    101 0.46 0.11 0.45 0.28  0.71 1687.61    1
    #>    102 0.48 0.16 0.46 0.24  0.84 1955.14    1
    #>    103 0.51 0.21 0.48 0.23  0.97 1776.23    1
    #>    104 0.54 0.24 0.50 0.24  1.07 1566.12    1
    #>    105 0.57 0.25 0.52 0.27  1.15 1757.57    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.15 0.76 0.51  1.10 1769.09    1
    #>    102 0.79 0.25 0.76 0.42  1.31 1816.33    1
    #>    103 0.87 0.36 0.82 0.43  1.68 1749.27    1
    #>    104 0.88 0.39 0.81 0.45  1.85 1742.88    1
    #>    105 0.92 0.46 0.82 0.47  2.16 2078.94    1
