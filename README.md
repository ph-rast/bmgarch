<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with rmarkdown::render("README.Rmd", output_format = "md_document") -->
<!-- badges: start -->
[![Travis build
status](https://travis-ci.com/ph-rast/bmgarch.svg?branch=master)](https://travis-ci.com/ph-rast/bmgarch)
![R-CMD-check](https://github.com/ph-rast/bmgarch/workflows/R-CMD-check/badge.svg)
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

`bmgarch` is not yet available on CRAN.

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
    #> Date:  Tue Sep  1 14:35:55 2020
    #> Elapsed time (min):  17.6
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> R_Ng-Ps 0.04 0.44 0.09 -0.89  0.89 1070.54    1
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.81 0.78 0.56 0.02  2.92  23.84 1.08
    #> var_Ng 1.28 0.35 1.37 0.38  1.92 176.82 1.03
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> A_Ps-Ps 0.35 0.09 0.36  0.15  0.53 1046.96 1.01
    #> A_Ng-Ps 0.04 0.08 0.03 -0.08  0.20    8.36 1.16
    #> A_Ps-Ng 0.06 0.12 0.10 -0.23  0.30  953.32 1.02
    #> A_Ng-Ng 0.46 0.15 0.45  0.16  0.66    3.53 1.50
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.70 0.21  0.77  0.12  0.94  65.07 1.07
    #> B_Ng-Ps -0.06 0.14 -0.04 -0.39  0.30 346.98 1.03
    #> B_Ps-Ng  0.34 0.38  0.38 -0.66  1.02  30.59 1.09
    #> B_Ng-Ng  0.30 0.19  0.24  0.03  0.72  11.71 1.13
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos  0.03 0.17  0.04 -0.35  0.29   7.93 1.18
    #> (Intercept)_Neg  0.11 0.11  0.13 -0.12  0.33  14.19 1.10
    #> Phi_Pos-Pos      0.08 0.38  0.11 -0.74  0.61   5.62 1.27
    #> Phi_Pos-Neg     -0.24 0.49 -0.35 -0.90  0.77   7.89 1.22
    #> Phi_Neg-Pos     -0.18 0.31 -0.27 -0.80  0.55 350.52 1.02
    #> Phi_Neg-Neg      0.09 0.39  0.22 -0.73  0.81 248.27 1.03
    #> Theta_Pos-Pos   -0.17 0.40 -0.21 -0.73  0.68   5.50 1.28
    #> Theta_Pos-Neg    0.18 0.54  0.25 -0.89  0.82   5.57 1.31
    #> Theta_Neg-Pos    0.23 0.32  0.32 -0.54  0.83  64.86 1.06
    #> Theta_Neg-Neg   -0.13 0.42 -0.27 -0.85  0.74  23.04 1.08
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   49.31   25.53   42.75   13.84  112.55 1101.46    1.01 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -800.77    5.75 -800.37 -813.18 -791.70    4.56    1.37

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.83 2.87 -0.85 -6.40  4.73  848.12 1.01
    #>    202 -0.53 2.75 -0.60 -5.63  4.70  613.11 1.02
    #>    203 -0.33 2.61 -0.31 -5.49  4.79  919.44 1.01
    #>    204 -0.31 2.57 -0.26 -5.39  4.79 1276.71 1.01
    #>    205 -0.13 2.58 -0.12 -5.27  4.91 1778.57 1.00
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.65 1.55 0.67 -2.44  3.60   57.49 1.04
    #>    202 0.46 1.62 0.43 -2.79  3.77  522.73 1.01
    #>    203 0.32 1.65 0.36 -2.89  3.49  728.16 1.01
    #>    204 0.31 1.65 0.26 -2.91  3.44 1805.40 1.00
    #>    205 0.25 1.57 0.25 -2.92  3.31 1546.44 1.00
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 7.32 2.59 6.46 3.99 13.49   36.97 1.06
    #>    202 6.47 3.40 5.71 3.15 14.41   28.52 1.07
    #>    203 6.23 4.40 5.16 2.98 14.78  101.40 1.03
    #>    204 6.08 4.33 4.94 2.93 15.13  443.08 1.02
    #>    205 6.06 4.48 4.81 2.87 15.42 1085.44 1.01
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.92 0.37 1.86 1.36  2.74 1108.64    1
    #>    202 2.24 1.00 2.00 1.41  4.58 1677.01    1
    #>    203 2.32 1.17 2.01 1.42  5.31 1992.77    1
    #>    204 2.38 1.61 2.01 1.42  5.66 1906.44    1
    #>    205 2.37 1.75 1.99 1.40  5.29 1951.10    1
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #>    201 -0.11 0.17 -0.12 -0.39  0.23  42.66 1.04
    #>    202 -0.05 0.21 -0.05 -0.47  0.38 368.15 1.02
    #>    203 -0.02 0.22 -0.01 -0.46  0.44  41.02 1.06
    #>    204  0.01 0.22  0.01 -0.44  0.49  31.00 1.07
    #>    205  0.01 0.22  0.02 -0.43  0.47  37.05 1.06

    plot(fit.fc, askNewPage = FALSE, type = "var")

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />


    plot(fit.fc, askNewPage = FALSE, type = "cor")

<img src="man/figures/README-forecastPlot-3.png" width="100%" />

Example 2: Stocks
-----------------

Here we use the first 100 days of Stata's stocks data on daily lagged
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
    #> Date:  Tue Sep  1 14:37:10 2020
    #> Elapsed time (min):  0.96
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.10 0.09 0.08 0.00  0.36 2378.97    1
    #> a_h_1,ns   0.08 0.07 0.07 0.00  0.27 2259.80    1
    #> a_h_1,hn   0.11 0.09 0.09 0.00  0.33 2183.70    1
    #> b_h_1,ty   0.45 0.18 0.47 0.10  0.77 1466.95    1
    #> b_h_1,ns   0.38 0.20 0.36 0.06  0.77 1176.61    1
    #> b_h_1,hn   0.39 0.17 0.38 0.08  0.73 1288.51    1
    #> c_h_var_ty 0.28 0.12 0.27 0.09  0.55 1283.87    1
    #> c_h_var_ns 0.35 0.13 0.35 0.12  0.63 1136.18    1
    #> c_h_var_hn 0.44 0.16 0.43 0.17  0.77 1250.73    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.65 0.06 0.65 0.51  0.75 1625.27    1
    #> R_hn-ty 0.73 0.05 0.74 0.62  0.82 2032.82    1
    #> R_hn-ns 0.64 0.06 0.64 0.50  0.75 2287.38    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.24  0.07 1351.67    1
    #> (Intercept)_nissan -0.01 0.08  0.00 -0.17  0.15 1496.11    1
    #> (Intercept)_honda  -0.02 0.09 -0.02 -0.19  0.16 1361.46    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   33.03   24.97   25.52    7.05   99.36 1809.21    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.50    5.17 -178.15 -189.69 -169.34  728.33    1.01

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
    #>    101 0.53 0.11 0.53 0.33  0.77 1805.69    1
    #>    102 0.58 0.16 0.56 0.34  0.98 1703.07    1
    #>    103 0.61 0.26 0.58 0.35  1.05 1680.22    1
    #>    104 0.63 0.25 0.59 0.36  1.17 1780.38    1
    #>    105 0.65 0.27 0.59 0.37  1.26 2153.90    1
    #>    106 0.65 0.25 0.60 0.37  1.21 1961.94    1
    #>    107 0.65 0.28 0.60 0.37  1.29 2052.69    1
    #>    108 0.66 0.30 0.60 0.38  1.30 2002.68    1
    #>    109 0.66 0.30 0.60 0.38  1.34 2007.72    1
    #>    110 0.66 0.30 0.60 0.38  1.41 2101.81    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.61 0.10 0.60 0.42  0.84 1660.92    1
    #>    102 0.64 0.17 0.61 0.43  0.98 1615.46    1
    #>    103 0.66 0.20 0.63 0.43  1.07 1725.27    1
    #>    104 0.67 0.22 0.63 0.44  1.10 1940.99    1
    #>    105 0.68 0.69 0.63 0.43  1.12 2000.14    1
    #>    106 0.73 2.68 0.63 0.43  1.15 2000.16    1
    #>    107 0.69 0.73 0.63 0.43  1.17 1949.52    1
    #>    108 0.68 0.31 0.63 0.43  1.20 1681.40    1
    #>    109 0.68 0.33 0.64 0.43  1.16 1823.36    1
    #>    110 0.68 0.45 0.63 0.43  1.18 1900.79    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.76 0.14 0.75 0.53  1.07 1567.26    1
    #>    102 0.83 0.23 0.79 0.53  1.44 1778.68    1
    #>    103 0.87 0.35 0.82 0.53  1.59 2039.45    1
    #>    104 0.89 0.32 0.83 0.55  1.72 2027.82    1
    #>    105 0.90 0.39 0.83 0.54  1.72 1953.13    1
    #>    106 0.91 0.44 0.83 0.55  1.72 1979.25    1
    #>    107 0.91 0.41 0.84 0.56  1.70 1887.18    1
    #>    108 0.91 0.39 0.83 0.55  1.78 1896.24    1
    #>    109 0.91 0.42 0.83 0.55  1.83 2094.37    1
    #>    110 0.92 0.42 0.84 0.55  1.71 1968.70    1

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
    #> Date:  Tue Sep  1 14:55:03 2020
    #> Elapsed time (min):  15.97
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.04 0.03  0.36  5.14 1.31
    #> a_h_1,ns   0.07 0.05 0.06 0.01  0.21  4.58 1.35
    #> a_h_1,hn   0.14 0.15 0.04 0.01  0.36  2.33 2.67
    #> b_h_1,ty   0.50 0.11 0.51 0.22  0.65  4.94 1.35
    #> b_h_1,ns   0.35 0.17 0.33 0.15  0.72  3.37 1.59
    #> b_h_1,hn   0.44 0.23 0.49 0.11  0.71  2.39 2.48
    #> c_h_var_ty 0.27 0.07 0.26 0.14  0.45  7.42 1.21
    #> c_h_var_ns 0.35 0.11 0.40 0.15  0.53  3.18 1.69
    #> c_h_var_hn 0.34 0.10 0.36 0.18  0.56  5.01 1.32
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> a_q 0.20 0.07 0.22 0.09  0.34  4.13 1.43
    #> b_q 0.16 0.12 0.12 0.04  0.43  3.27 1.63
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> S_ns-ty 0.63 0.05 0.63 0.49  0.73 13.47 1.17
    #> S_hn-ty 0.71 0.05 0.69 0.63  0.81  4.60 1.35
    #> S_hn-ns 0.67 0.06 0.70 0.52  0.74  4.43 1.38
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #> (Intercept)_toyota  -0.02 0.14 -0.02 -0.20  0.18  2.17 3.54
    #> (Intercept)_nissan   0.13 0.08  0.18 -0.09  0.20  2.62 2.04
    #> (Intercept)_honda    0.02 0.12  0.03 -0.14  0.17  2.59 2.16
    #> Phi_toyota-toyota    0.04 0.42  0.20 -0.58  0.54  2.26 2.97
    #> Phi_toyota-nissan   -0.16 0.51 -0.52 -0.65  0.56  2.29 2.81
    #> Phi_toyota-honda    -0.15 0.53 -0.40 -0.74  0.58  2.25 3.03
    #> Phi_nissan-toyota    0.36 0.35  0.50 -0.33  0.83  3.04 1.81
    #> Phi_nissan-nissan   -0.37 0.26 -0.45 -0.62  0.37  4.20 1.42
    #> Phi_nissan-honda    -0.06 0.32 -0.08 -0.53  0.60  3.15 1.70
    #> Phi_honda-toyota    -0.24 0.27 -0.36 -0.69  0.33  4.13 1.43
    #> Phi_honda-nissan     0.12 0.29  0.25 -0.38  0.70  3.87 1.48
    #> Phi_honda-honda     -0.34 0.37 -0.56 -0.70  0.29  2.47 2.31
    #> Theta_toyota-toyota -0.19 0.45 -0.21 -0.82  0.37  2.31 2.77
    #> Theta_toyota-nissan  0.27 0.50  0.63 -0.42  0.81  2.30 2.81
    #> Theta_toyota-honda   0.26 0.63  0.70 -0.56  0.89  2.16 3.73
    #> Theta_nissan-toyota -0.43 0.49 -0.75 -0.91  0.30  2.41 2.51
    #> Theta_nissan-nissan  0.39 0.23  0.43 -0.34  0.59  5.22 1.32
    #> Theta_nissan-honda   0.07 0.37  0.03 -0.63  0.65  2.71 2.00
    #> Theta_honda-toyota  -0.04 0.41  0.03 -0.64  0.47  2.52 2.21
    #> Theta_honda-nissan  -0.02 0.31  0.13 -0.65  0.54  4.00 1.46
    #> Theta_honda-honda    0.52 0.45  0.87 -0.21  0.99  2.36 2.60
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  52.04  20.12  48.68  14.04 100.98   5.67   1.30 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -170.53    5.19 -170.09 -183.93 -161.78    3.11    1.76
    fc <- forecast(fit3, ahead =  10)

    plot( fc,askNewPage = FALSE, type =  'mean' ) 

<img src="man/figures/README-fit3ForecastPlot-1.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-2.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-3.png" width="100%" />

### Compute Model Weights

Obtain model weights with either the stacking or the pseudo BMA method.
These methods are inherited from the `loo` package.

First, gather models to a `bmgarch_list`.

    ## use bmgarch_list function to collect bmgarch objects
    modfits <- suppressMessages( bmgarch_list(fit1, fit2, fit3) )

Compute model weights with the stacking method (default) and the the
approximate (default) leave-future-out cross validation (LFO CV). `L`
defines the minimal length of the time series before we start engaging
in cross-validation. Eg., for a time series with length 100, `L = 50`
reserves values 51--100 as the cross-validation sample. Note that the
standard is to use the approximate `backward` method to CV as it results
in fewest refits. Exact CV is also available with `exact` but not
encouraged as it results in refitting all CV models.

    mw <- suppressMessages( model_weights(modfits, L = 80, method = 'stacking' ) )
    #> 
    #> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
    #> 
    #> COMPILING MODEL 'CCCMGARCH' NOW.
    #> 
    #> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
    #> Using threshold  0.6 , model was refit  1  times, at observations 80 
    #> Using threshold  0.6 , model was refit  0  times, at observations 
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
    #> Using threshold  0.6 , model was refit  3  times, at observations 96 87 84

    ## Return model weights:
    mw
    #> Method: stacking
    #> ------
    #>        weight
    #> model1 0.000 
    #> model2 1.000 
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
    #>    101 -0.09 0.75 -0.08 -1.54  1.39    NA   NA
    #>    102 -0.09 0.78 -0.11 -1.59  1.42    NA   NA
    #>    103 -0.08 0.80 -0.09 -1.66  1.44    NA   NA
    #>    104 -0.09 0.80 -0.09 -1.69  1.48    NA   NA
    #>    105 -0.11 0.82 -0.10 -1.81  1.51    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101  0.03 0.84  0.03 -1.60  1.65    NA   NA
    #>    102  0.03 0.87  0.00 -1.64  1.83    NA   NA
    #>    103 -0.01 0.88  0.01 -1.81  1.70    NA   NA
    #>    104  0.00 0.89 -0.01 -1.80  1.71    NA   NA
    #>    105 -0.02 0.89 -0.03 -1.75  1.70    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.01 0.93 -0.02 -1.89  1.85    NA   NA
    #>    102 -0.04 0.98 -0.05 -1.99  1.87    NA   NA
    #>    103  0.00 1.00  0.01 -1.97  1.95    NA   NA
    #>    104 -0.05 1.00 -0.05 -2.03  2.00    NA   NA
    #>    105 -0.03 0.98  0.00 -1.99  1.84    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.52 0.11 0.51 0.33  0.76    NA   NA
    #>    102 0.54 0.16 0.53 0.33  0.82    NA   NA
    #>    103 0.58 0.17 0.55 0.34  0.97    NA   NA
    #>    104 0.60 0.20 0.57 0.34  1.03    NA   NA
    #>    105 0.61 0.21 0.58 0.35  1.14    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.63 0.12 0.62 0.44  0.90    NA   NA
    #>    102 0.64 0.15 0.62 0.43  0.98    NA   NA
    #>    103 0.67 0.24 0.63 0.43  1.11    NA   NA
    #>    104 0.68 0.25 0.64 0.43  1.15    NA   NA
    #>    105 0.68 0.22 0.64 0.43  1.18    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.77 0.14 0.76 0.54  1.07    NA   NA
    #>    102 0.78 0.18 0.76 0.50  1.19    NA   NA
    #>    103 0.86 0.40 0.80 0.52  1.56    NA   NA
    #>    104 0.88 0.40 0.82 0.51  1.66    NA   NA
    #>    105 0.92 0.48 0.83 0.52  1.85    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.06 0.65  0.5  0.75    NA   NA
    #>    102 0.64 0.06 0.65  0.5  0.75    NA   NA
    #>    103 0.64 0.06 0.65  0.5  0.75    NA   NA
    #>    104 0.64 0.06 0.65  0.5  0.75    NA   NA
    #>    105 0.64 0.06 0.65  0.5  0.75    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    102 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    103 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    104 0.73 0.05 0.73 0.62  0.82    NA   NA
    #>    105 0.73 0.05 0.73 0.62  0.82    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.07 0.64  0.5  0.76    NA   NA
    #>    102 0.64 0.07 0.64  0.5  0.76    NA   NA
    #>    103 0.64 0.07 0.64  0.5  0.76    NA   NA
    #>    104 0.64 0.07 0.64  0.5  0.76    NA   NA
    #>    105 0.64 0.07 0.64  0.5  0.76    NA   NA

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
    #> Date:  Tue Sep  1 15:39:48 2020
    #> Elapsed time (min):  0.66
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.06 0.00  0.33 1652.43    1
    #> a_h_1,hn   0.08 0.08 0.06 0.00  0.30 2214.30    1
    #> a_h_2,ty   0.09 0.09 0.07 0.00  0.31 2200.42    1
    #> a_h_2,hn   0.12 0.12 0.08 0.00  0.44 1930.74    1
    #> b_h_1,ty   0.21 0.17 0.17 0.01  0.60 2010.68    1
    #> b_h_1,hn   0.18 0.15 0.15 0.01  0.53 1577.43    1
    #> b_h_2,ty   0.26 0.17 0.24 0.01  0.62 2062.16    1
    #> b_h_2,hn   0.20 0.17 0.15 0.01  0.63 1770.49    1
    #> c_h_var_ty 0.22 0.10 0.21 0.07  0.47 1458.45    1
    #> c_h_var_hn 0.39 0.15 0.38 0.14  0.72 1254.30    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.81 2553.57    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.26  0.07 1749.19    1
    #> (Intercept)_honda  -0.05 0.09 -0.05 -0.23  0.13 1586.00    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.61 0.49 -1.57 -2.65 -0.75 1404.15    1
    #> beta0_hn -1.03 0.44 -0.97 -1.99 -0.32 1115.48    1
    #> beta_ty  -0.18 0.39 -0.18 -0.94  0.57 1485.92    1
    #> beta_hn   0.04 0.33  0.06 -0.69  0.65 1319.36    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   45.42   28.53   38.77    9.41  116.20 3323.47    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -129.85    4.14 -129.53 -138.57 -122.82  755.26    1.00

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
    #>    101 0.46 0.11 0.46 0.28  0.70 1941.21    1
    #>    102 0.48 0.16 0.46 0.24  0.82 2163.48    1
    #>    103 0.52 0.22 0.48 0.23  1.03 2177.91    1
    #>    104 0.54 0.24 0.50 0.24  1.08 2249.09    1
    #>    105 0.57 0.26 0.53 0.27  1.15 1949.69    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.14 0.76 0.52  1.09 2072.75    1
    #>    102 0.79 0.24 0.76 0.44  1.29 2241.84    1
    #>    103 0.87 0.34 0.81 0.44  1.73 2180.81    1
    #>    104 0.90 0.50 0.81 0.46  1.98 2069.16    1
    #>    105 0.93 0.59 0.82 0.49  2.14 2006.15    1
