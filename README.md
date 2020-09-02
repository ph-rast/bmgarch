<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with rmarkdown::render("README.Rmd", output_format = "md_document") -->
<!-- badges: start -->
[![Travis build
status](https://travis-ci.com/ph-rast/bmgarch.svg?branch=master)](https://travis-ci.com/ph-rast/bmgarch)
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
    #> Date:  Tue Sep  1 22:06:52 2020
    #> Elapsed time (min):  14.78
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean   sd mdn  2.5% 97.5%   n_eff Rhat
    #> R_Ng-Ps 0.01 0.49   0 -0.88   0.9 1225.16    1
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.77 0.98 0.34 0.01  3.53 277.21    1
    #> var_Ng 1.20 0.45 1.27 0.16  1.95 361.59    1
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #> A_Ps-Ps 0.34 0.10 0.34  0.13  0.55 1236.60    1
    #> A_Ng-Ps 0.06 0.08 0.06 -0.09  0.21 1539.86    1
    #> A_Ps-Ng 0.05 0.14 0.06 -0.23  0.34 1001.81    1
    #> A_Ng-Ng 0.39 0.12 0.39  0.13  0.61  984.89    1
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.71 0.25  0.83  0.05  0.95 232.08 1.01
    #> B_Ng-Ps -0.09 0.18 -0.09 -0.52  0.33 252.72 1.01
    #> B_Ps-Ng  0.27 0.40  0.27 -0.67  1.16 388.01 1.01
    #> B_Ng-Ng  0.35 0.21  0.34  0.02  0.78 785.91 1.00
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos -0.03 0.17 -0.02 -0.42  0.29 683.62 1.00
    #> (Intercept)_Neg  0.07 0.12  0.06 -0.16  0.34 949.69 1.01
    #> Phi_Pos-Pos     -0.06 0.39 -0.03 -0.87  0.67 577.24 1.00
    #> Phi_Pos-Neg     -0.10 0.46 -0.12 -0.89  0.81 534.97 1.00
    #> Phi_Neg-Pos     -0.16 0.36 -0.16 -0.85  0.62 626.77 1.01
    #> Phi_Neg-Neg      0.05 0.43  0.06 -0.71  0.82 858.42 1.00
    #> Theta_Pos-Pos   -0.02 0.41 -0.04 -0.79  0.82 574.07 1.00
    #> Theta_Pos-Neg    0.01 0.47  0.05 -0.92  0.79 506.37 1.00
    #> Theta_Neg-Pos    0.18 0.36  0.18 -0.61  0.86 595.10 1.02
    #> Theta_Neg-Neg   -0.06 0.44 -0.06 -0.86  0.74 810.21 1.00
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   52.18   28.07   47.06   13.69  119.03 2529.13    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -803.41    4.97 -803.11 -814.06 -794.68  383.58    1.00

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.25 2.91 -0.29 -6.00  5.46 1860.08    1
    #>    202 -0.19 2.77 -0.20 -5.60  5.41 2037.99    1
    #>    203 -0.13 2.75 -0.08 -5.68  4.97 1952.89    1
    #>    204 -0.14 2.59 -0.16 -5.24  5.03 1952.77    1
    #>    205 -0.09 2.60 -0.09 -5.20  5.28 1870.26    1
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.35 1.45 0.36 -2.53  3.16 1698.74    1
    #>    202 0.25 1.58 0.22 -2.72  3.54 1594.92    1
    #>    203 0.15 1.59 0.14 -3.14  3.26 1923.31    1
    #>    204 0.20 1.61 0.25 -3.01  3.33 2060.69    1
    #>    205 0.14 1.57 0.10 -2.98  3.24 2014.93    1
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 7.77 2.61 7.48 3.89 13.49  570.52    1
    #>    202 7.16 3.28 6.52 3.36 15.00  898.08    1
    #>    203 6.73 3.73 5.82 3.11 15.35 1172.30    1
    #>    204 6.51 4.63 5.33 3.04 16.12 1578.67    1
    #>    205 6.29 4.65 5.08 2.93 16.37 1613.58    1
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.97 0.43 1.90 1.36  3.02 1031.77    1
    #>    202 2.22 0.76 2.03 1.40  4.11 2006.64    1
    #>    203 2.31 0.99 2.07 1.41  4.62 2062.82    1
    #>    204 2.34 1.00 2.08 1.41  4.97 2096.10    1
    #>    205 2.31 1.08 2.05 1.40  4.72 2017.50    1
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.07 0.17 -0.06 -0.40  0.23 1012.83    1
    #>    202 -0.07 0.21 -0.07 -0.46  0.37 1326.89    1
    #>    203 -0.06 0.21 -0.06 -0.50  0.39 1619.66    1
    #>    204 -0.04 0.21 -0.04 -0.46  0.41 1721.83    1
    #>    205 -0.03 0.20 -0.04 -0.45  0.41 1537.11    1

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
    #> Date:  Tue Sep  1 22:08:07 2020
    #> Elapsed time (min):  0.95
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.11 0.10 0.08 0.00  0.37 2185.86    1
    #> a_h_1,ns   0.08 0.07 0.06 0.00  0.26 2744.57    1
    #> a_h_1,hn   0.11 0.09 0.09 0.00  0.33 2377.30    1
    #> b_h_1,ty   0.45 0.18 0.46 0.10  0.75 1616.31    1
    #> b_h_1,ns   0.37 0.20 0.35 0.06  0.77 1373.82    1
    #> b_h_1,hn   0.39 0.18 0.38 0.08  0.76 1316.46    1
    #> c_h_var_ty 0.29 0.12 0.27 0.10  0.56 1455.06    1
    #> c_h_var_ns 0.35 0.13 0.36 0.10  0.60 1527.56    1
    #> c_h_var_hn 0.44 0.16 0.43 0.16  0.78 1264.84    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_ns-ty 0.65 0.07 0.65 0.50  0.76 2426.85    1
    #> R_hn-ty 0.73 0.05 0.74 0.62  0.82 2109.22    1
    #> R_hn-ns 0.64 0.07 0.64 0.50  0.75 2336.69    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.24  0.07 1729.36    1
    #> (Intercept)_nissan -0.01 0.08 -0.01 -0.16  0.15 1693.25    1
    #> (Intercept)_honda  -0.02 0.09 -0.02 -0.20  0.16 1574.42    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   31.68   23.48   24.19    7.09   91.31 2270.38    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.37    5.22 -178.00 -189.56 -169.32  705.81    1.00

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
    #>    101 0.54 0.11 0.53 0.34  0.77 1652.55    1
    #>    102 0.58 0.17 0.56 0.35  0.96 1983.69    1
    #>    103 0.61 0.22 0.58 0.36  1.09 1950.83    1
    #>    104 0.63 0.24 0.59 0.37  1.15 1862.93    1
    #>    105 0.65 0.29 0.60 0.37  1.27 2011.49    1
    #>    106 0.66 0.36 0.60 0.37  1.38 1816.63    1
    #>    107 0.67 0.33 0.61 0.38  1.38 1910.15    1
    #>    108 0.67 0.37 0.61 0.38  1.36 1816.65    1
    #>    109 0.68 0.32 0.61 0.38  1.35 2084.51    1
    #>    110 0.67 0.33 0.61 0.39  1.31 1920.16    1
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.61 0.10 0.60 0.43  0.84 2066.22    1
    #>    102 0.64 0.15 0.62 0.43  0.97 1752.66    1
    #>    103 0.66 0.18 0.62 0.43  1.05 1942.56    1
    #>    104 0.67 0.21 0.63 0.43  1.11 1810.40    1
    #>    105 0.67 0.23 0.63 0.43  1.20 1880.72    1
    #>    106 0.67 0.24 0.63 0.43  1.22 1737.02    1
    #>    107 0.68 0.23 0.63 0.43  1.23 1954.45    1
    #>    108 0.68 0.23 0.63 0.43  1.23 1931.12    1
    #>    109 0.68 0.27 0.63 0.44  1.16 1635.07    1
    #>    110 0.68 0.23 0.64 0.43  1.19 2138.37    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.76 0.14 0.75 0.53  1.07 1787.36    1
    #>    102 0.82 0.21 0.79 0.53  1.29 1905.72    1
    #>    103 0.86 0.26 0.81 0.54  1.44 1636.53    1
    #>    104 0.89 0.34 0.82 0.55  1.58 1658.82    1
    #>    105 0.92 0.70 0.83 0.55  1.67 1346.66    1
    #>    106 0.93 1.16 0.84 0.56  1.78 1752.68    1
    #>    107 0.96 2.27 0.84 0.56  1.79 1920.54    1
    #>    108 0.95 1.43 0.84 0.56  1.84 1973.01    1
    #>    109 0.94 0.59 0.84 0.55  1.85 2077.14    1
    #>    110 0.95 0.88 0.84 0.56  1.98 2018.69    1

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
    #> Date:  Tue Sep  1 22:24:48 2020
    #> Elapsed time (min):  15.41
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_h_1,ty   0.15 0.13 0.14 0.02  0.44   4.90 1.31
    #> a_h_1,ns   0.07 0.07 0.07 0.01  0.27   8.17 1.17
    #> a_h_1,hn   0.16 0.11 0.11 0.01  0.33   3.41 1.54
    #> b_h_1,ty   0.42 0.12 0.42 0.15  0.70  49.76 1.05
    #> b_h_1,ns   0.54 0.19 0.62 0.09  0.80   4.75 1.32
    #> b_h_1,hn   0.48 0.14 0.46 0.15  0.79 155.25 1.04
    #> c_h_var_ty 0.33 0.12 0.29 0.11  0.50   3.50 1.52
    #> c_h_var_ns 0.26 0.11 0.25 0.10  0.54   6.86 1.21
    #> c_h_var_hn 0.32 0.13 0.29 0.14  0.69  10.34 1.13
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> a_q 0.29 0.15 0.22 0.05  0.51  2.49 2.20
    #> b_q 0.29 0.11 0.30 0.04  0.48  7.21 1.18
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #> S_ns-ty 0.64 0.08 0.62 0.46  0.74  3.52 1.51
    #> S_hn-ty 0.72 0.07 0.74 0.61  0.82  3.56 1.50
    #> S_hn-ns 0.59 0.08 0.61 0.49  0.74  4.04 1.40
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota  -0.06 0.07 -0.07 -0.22  0.07  12.86 1.12
    #> (Intercept)_nissan   0.11 0.12  0.11 -0.13  0.28   2.82 1.84
    #> (Intercept)_honda   -0.02 0.09 -0.01 -0.19  0.21 116.92 1.04
    #> Phi_toyota-toyota   -0.01 0.36  0.04 -0.63  0.63   3.92 1.44
    #> Phi_toyota-nissan   -0.17 0.40 -0.02 -0.69  0.62   3.66 1.49
    #> Phi_toyota-honda     0.23 0.32  0.20 -0.56  0.70   6.65 1.23
    #> Phi_nissan-toyota    0.47 0.33  0.61 -0.44  0.87   9.28 1.17
    #> Phi_nissan-nissan   -0.32 0.29 -0.46 -0.76  0.45  11.13 1.15
    #> Phi_nissan-honda    -0.02 0.41  0.10 -0.63  0.72   3.93 1.45
    #> Phi_honda-toyota    -0.17 0.38 -0.25 -0.90  0.45   4.93 1.32
    #> Phi_honda-nissan    -0.06 0.42 -0.01 -0.57  0.86   3.77 1.46
    #> Phi_honda-honda     -0.29 0.31 -0.37 -0.71  0.47   5.58 1.28
    #> Theta_toyota-toyota -0.16 0.37 -0.14 -0.76  0.52   4.45 1.37
    #> Theta_toyota-nissan  0.34 0.41  0.19 -0.50  0.87   3.51 1.52
    #> Theta_toyota-honda  -0.14 0.33 -0.08 -0.68  0.56   5.21 1.30
    #> Theta_nissan-toyota -0.49 0.34 -0.64 -0.88  0.47   7.64 1.20
    #> Theta_nissan-nissan  0.34 0.28  0.46 -0.37  0.74   8.30 1.18
    #> Theta_nissan-honda  -0.02 0.43 -0.15 -0.80  0.56   3.34 1.59
    #> Theta_honda-toyota  -0.20 0.40 -0.06 -0.77  0.64   4.32 1.38
    #> Theta_honda-nissan   0.13 0.41  0.06 -0.79  0.69   4.71 1.34
    #> Theta_honda-honda    0.57 0.44  0.89 -0.45  0.91   3.28 1.61
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  36.72  22.41  35.72  11.66 100.92   8.59   1.15 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -171.51    8.38 -172.71 -187.27 -158.23    2.50    2.24
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
reserves values 51--100 as the cross-validation sample. Note that the
standard is to use the approximate `backward` method to CV as it results
in fewest refits. Exact CV is also available with `exact` but not
encouraged as it results in refitting all CV models.

    mw <- model_weights(modfits, L = 50, method = 'stacking')
    #> Using threshold  0.6 , model was refit  10  times, at observations 93 85 83 75 74 72 64 60 58 52

    ## Return model weights:
    mw
    #> Method: stacking
    #> ------
    #>        weight
    #> model1 0.555 
    #> model2 0.445 
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
    #>    101 -0.09 0.57 -0.09 -1.23  1.01    NA   NA
    #>    102 -0.08 0.59 -0.09 -1.23  1.13    NA   NA
    #>    103 -0.10 0.62 -0.09 -1.35  1.13    NA   NA
    #>    104 -0.10 0.61 -0.09 -1.37  1.11    NA   NA
    #>    105 -0.08 0.62 -0.07 -1.31  1.13    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101  0.00 0.62 -0.01 -1.22  1.22    NA   NA
    #>    102  0.00 0.65  0.02 -1.31  1.24    NA   NA
    #>    103 -0.01 0.64  0.00 -1.26  1.26    NA   NA
    #>    104  0.00 0.64 -0.01 -1.33  1.29    NA   NA
    #>    105  0.01 0.63  0.02 -1.21  1.27    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.04 0.69 -0.04 -1.39  1.31    NA   NA
    #>    102 -0.02 0.70  0.00 -1.39  1.35    NA   NA
    #>    103 -0.05 0.74 -0.03 -1.54  1.40    NA   NA
    #>    104 -0.02 0.75 -0.02 -1.53  1.50    NA   NA
    #>    105  0.00 0.73  0.00 -1.44  1.43    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.53 0.08 0.53 0.39  0.70    NA   NA
    #>    102 0.57 0.11 0.55 0.40  0.82    NA   NA
    #>    103 0.60 0.16 0.58 0.41  0.97    NA   NA
    #>    104 0.62 0.18 0.59 0.41  1.04    NA   NA
    #>    105 0.64 0.20 0.60 0.42  1.13    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.62 0.08 0.61 0.48  0.80    NA   NA
    #>    102 0.64 0.12 0.63 0.48  0.89    NA   NA
    #>    103 0.66 0.14 0.64 0.48  1.00    NA   NA
    #>    104 0.67 0.16 0.64 0.48  1.02    NA   NA
    #>    105 0.67 0.16 0.65 0.48  1.04    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.77 0.10 0.76 0.60  1.00    NA   NA
    #>    102 0.81 0.16 0.79 0.58  1.14    NA   NA
    #>    103 0.87 0.22 0.83 0.61  1.38    NA   NA
    #>    104 0.89 0.25 0.84 0.61  1.50    NA   NA
    #>    105 0.91 0.28 0.86 0.62  1.57    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.05 0.65 0.55  0.73    NA   NA
    #>    102 0.65 0.05 0.65 0.55  0.73    NA   NA
    #>    103 0.65 0.05 0.65 0.55  0.73    NA   NA
    #>    104 0.65 0.05 0.65 0.55  0.73    NA   NA
    #>    105 0.65 0.05 0.65 0.55  0.73    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.04 0.74 0.65   0.8    NA   NA
    #>    102 0.73 0.04 0.74 0.65   0.8    NA   NA
    #>    103 0.73 0.04 0.74 0.65   0.8    NA   NA
    #>    104 0.73 0.04 0.74 0.65   0.8    NA   NA
    #>    105 0.73 0.04 0.74 0.65   0.8    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.05 0.64 0.54  0.72    NA   NA
    #>    102 0.64 0.05 0.64 0.54  0.72    NA   NA
    #>    103 0.64 0.05 0.64 0.54  0.72    NA   NA
    #>    104 0.64 0.05 0.64 0.54  0.72    NA   NA
    #>    105 0.64 0.05 0.64 0.54  0.72    NA   NA

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
    #> Date:  Wed Sep  2 00:16:12 2020
    #> Elapsed time (min):  0.58
    #> 
    #> GARCH(2,2)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> a_h_1,ty   0.09 0.09 0.07 0.00  0.33 1786.48    1
    #> a_h_1,hn   0.08 0.08 0.05 0.00  0.29 1964.84    1
    #> a_h_2,ty   0.10 0.09 0.07 0.00  0.32 2107.94    1
    #> a_h_2,hn   0.13 0.13 0.09 0.00  0.50 1681.02    1
    #> b_h_1,ty   0.21 0.16 0.17 0.01  0.58 2000.09    1
    #> b_h_1,hn   0.18 0.14 0.14 0.00  0.54 2162.12    1
    #> b_h_2,ty   0.26 0.17 0.24 0.01  0.62 1974.25    1
    #> b_h_2,hn   0.18 0.16 0.15 0.01  0.58 1567.39    1
    #> c_h_var_ty 0.22 0.10 0.20 0.07  0.44 1408.81    1
    #> c_h_var_hn 0.40 0.15 0.39 0.14  0.73 1542.60    1
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #> R_hn-ty 0.73 0.05 0.73 0.62  0.81 2491.25    1
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> (Intercept)_toyota -0.09 0.08 -0.09 -0.26  0.07 1915.77    1
    #> (Intercept)_honda  -0.05 0.10 -0.05 -0.25  0.13 1931.59    1
    #> 
    #> 
    #> Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):
    #> 
    #>           mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #> beta0_ty -1.63 0.48 -1.60 -2.64 -0.81 1272.66    1
    #> beta0_hn -1.00 0.43 -0.94 -1.95 -0.32 1244.65    1
    #> beta_ty  -0.18 0.38 -0.19 -0.92  0.61 2001.59    1
    #> beta_hn   0.06 0.32  0.08 -0.62  0.65 1630.81    1
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #>   45.92   28.82   39.54    8.94  116.75 2758.71    1.00 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -129.67    4.07 -129.39 -138.46 -122.67  798.37    1.00

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
    #>    101 0.46 0.11 0.45 0.27  0.71 1890.91    1
    #>    102 0.47 0.16 0.45 0.24  0.82 1777.67    1
    #>    103 0.51 0.21 0.47 0.23  1.00 1986.70    1
    #>    104 0.53 0.24 0.50 0.23  1.11 1880.80    1
    #>    105 0.56 0.25 0.52 0.26  1.13 1941.90    1
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    101 0.77 0.14 0.76 0.52  1.09 2052.65    1
    #>    102 0.80 0.25 0.77 0.43  1.33 1949.45    1
    #>    103 0.89 0.39 0.83 0.43  1.81 1950.88    1
    #>    104 0.90 0.41 0.82 0.43  1.86 2028.12    1
    #>    105 0.94 0.76 0.83 0.46  2.08 1969.39    1
