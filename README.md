
<!-- README.md is generated from README.Rmd. Please edit that file -->
bmgarch
=======

`bmgarch` estimates Bayesian multivariate generalized autoregressive conditional heteroskedasticity (MGARCH) models. Currently, bmgarch supports ARMA(1,1) and intercept-only (Constant) mean structures, and a variety of MGARCH(P,Q) parameterizations. In increasing order of complexity:

-   CCC(P, Q): Constant Conditional Correlation
-   DCC(P, Q): Dynamic Conditional Correlation
-   BEKK(P, Q): Baba, Engle, Kraft, and Kroner
-   pdBEKK(P, Q): BEKK(P, Q) with positive diagonal constraints

Installation
------------

`bmgarch` is not yet available on CRAN.

The development version can be installed from [GitHub](https://github.com/) with:

``` r
devtools::install_github("ph-rast/bmgarch")
```

Example 1: Behavioral Data
--------------------------

In this example, we use the pdBEKK(1,1) model for the variances, and an intercept-only model for the means.

``` r
library(bmgarch)

data(panas)
head(panas)
```

|     Pos|     Neg|
|-------:|-------:|
|  -2.193|  -2.419|
|   1.567|  -0.360|
|  -0.124|  -1.202|
|   0.020|  -1.311|
|  -0.150|   2.004|
|   3.877|   1.008|

``` r

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
```

### Parameter estimates

``` r
summary(fit)
#> Model: pdBEKK-MGARCH
#> Basic Specification: H_t = D_t R D_t
#> H_t = C + A'[y_(t-1)*y'_(t-1)]A + B'H_(t-1)B
#> 
#> Distribution:  Student_t
#> ---
#> Iterations:  1000
#> Chains:  4
#> Date:  Fri Aug 21 17:48:50 2020
#> Elapsed time (min):  15.6
#> 
#> ---
#> Constant correlation, R (diag[C]*R*diag[C]):
#> 
#>          mean   sd  mdn  2.5% 97.5%  n_eff Rhat
#> R_Ng-Ps -0.03 0.45 -0.1 -0.88  0.89 957.21 1.01
#> 
#> 
#> Constant variances (diag[C]):
#> 
#>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> var_Ps 0.49 0.79 0.16  0.0  3.13  19.85 1.10
#> var_Ng 1.28 0.37 1.35  0.2  1.89 172.21 1.04
#> 
#> 
#> MGARCH(1,1) estimates for A:
#> 
#>         mean   sd  mdn  2.5% 97.5%  n_eff Rhat
#> A_Ps-Ps 0.37 0.11 0.38  0.15  0.53   6.70 1.21
#> A_Ng-Ps 0.05 0.07 0.03 -0.08  0.21  14.32 1.10
#> A_Ps-Ng 0.03 0.12 0.01 -0.20  0.29 321.43 1.03
#> A_Ng-Ng 0.39 0.11 0.38  0.15  0.61 770.80 1.01
#> 
#> 
#> MGARCH(1,1) estimates for B:
#> 
#>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#> B_Ps-Ps  0.75 0.21  0.82  0.10  0.94 153.93 1.04
#> B_Ng-Ps -0.07 0.15 -0.05 -0.40  0.27 112.76 1.04
#> B_Ps-Ng  0.31 0.36  0.38 -0.69  1.04 283.69 1.02
#> B_Ng-Ng  0.33 0.18  0.29  0.02  0.71 824.59 1.02
#> 
#> 
#> ARMA(1,1) estimates on the location:
#> 
#>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#> (Intercept)_Pos  0.04 0.17  0.04 -0.33  0.27   5.28 1.28
#> (Intercept)_Neg  0.06 0.11  0.03 -0.14  0.34 226.10 1.03
#> Phi_Pos-Pos      0.04 0.35  0.13 -0.76  0.65  19.15 1.10
#> Phi_Pos-Neg     -0.31 0.52 -0.34 -0.91  0.72   4.50 1.35
#> Phi_Neg-Pos     -0.10 0.32 -0.02 -0.78  0.58  64.70 1.05
#> Phi_Neg-Neg      0.24 0.53  0.25 -0.80  0.89   4.10 1.41
#> Theta_Pos-Pos   -0.15 0.39 -0.22 -0.75  0.69  10.18 1.15
#> Theta_Pos-Neg    0.22 0.53  0.26 -0.83  0.84   4.64 1.34
#> Theta_Neg-Pos    0.14 0.31  0.05 -0.55  0.80 274.54 1.03
#> Theta_Neg-Neg   -0.27 0.57 -0.26 -0.98  0.82   3.78 1.47
#> 
#> 
#> Df constant student_t (nu):
#> 
#>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
#>  62.00  29.08  61.23  15.46 111.27   5.99   1.22 
#> 
#> 
#> Log density posterior estimate:
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#> -801.25    5.62 -801.14 -813.07 -792.30    4.63    1.38
```

### Forecasted values

``` r
fit.fc <- forecast(fit, ahead = 5)

fit.fc
#> ---
#> [Mean] Forecast for 5 ahead:
#> 
#> Pos :
#>       
#> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
#>    201 -0.79 3.28 -0.77 -7.85  5.70  204.02 1.02
#>    202 -0.60 3.06 -0.61 -6.99  5.21  306.47 1.02
#>    203 -0.41 3.01 -0.37 -6.46  5.42  922.14 1.02
#>    204 -0.33 2.92 -0.26 -6.46  5.29 1818.66 1.01
#>    205 -0.14 2.84 -0.18 -5.77  5.63 1694.03 1.00
#> Neg :
#>       
#> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
#>    201 0.62 1.59 0.60 -2.54  3.66   47.31 1.04
#>    202 0.49 1.56 0.51 -2.89  3.38  184.78 1.02
#>    203 0.43 1.60 0.38 -2.81  3.63  418.74 1.01
#>    204 0.27 1.52 0.29 -2.61  3.20 1391.86 1.01
#>    205 0.23 1.55 0.25 -2.94  3.25 1775.12 1.00
#> ---
#> [Variance] Forecast for 5 ahead:
#> 
#> Pos :
#>       
#> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#>    201 9.27 3.51 9.00 4.10 16.99   8.38 1.17
#>    202 8.44 4.60 7.62 3.44 19.13  28.33 1.05
#>    203 7.86 4.91 6.65 3.21 20.29  43.53 1.04
#>    204 7.57 5.54 6.05 3.03 20.93 105.40 1.03
#>    205 7.27 5.71 5.63 2.93 20.77 197.20 1.02
#> Neg :
#>       
#> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#>    201 1.89 0.37 1.81 1.39  2.85  75.58 1.03
#>    202 2.13 0.80 1.95 1.42  3.89 466.55 1.02
#>    203 2.19 1.06 1.96 1.42  4.33 397.26 1.01
#>    204 2.18 0.88 1.94 1.42  4.18 745.91 1.01
#>    205 2.18 1.04 1.94 1.43  4.22 352.44 1.01
#> [Correlation] Forecast for 5 ahead:
#> 
#> Neg_Pos :
#>       
#> period  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#>    201 -0.09 0.15 -0.11 -0.37  0.22 157.14 1.02
#>    202 -0.06 0.19 -0.05 -0.43  0.31 827.46 1.00
#>    203 -0.03 0.20 -0.02 -0.41  0.38 545.83 1.01
#>    204 -0.02 0.20 -0.01 -0.41  0.38 428.46 1.01
#>    205 -0.01 0.20  0.00 -0.43  0.38 930.12 1.01
```

``` r
plot(fit.fc, askNewPage = FALSE, type = "var")
```

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />

``` r

plot(fit.fc, askNewPage = FALSE, type = "cor")
```

<img src="man/figures/README-forecastPlot-3.png" width="100%" />

Example 2: Stocks
-----------------

Here we use the first 100 days of Stata's stocks data on daily lagged returns of three Japanese automakers, Toyota, Nissan, and Honda.

``` r
library(bmgarch)

data(stocks)
head(stocks)
```

| date       |    t|      toyota|      nissan|       honda|
|:-----------|----:|-----------:|-----------:|-----------:|
| 2003-01-02 |    1|   0.0151675|   0.0294704|   0.0316103|
| 2003-01-03 |    2|   0.0048201|   0.0081735|   0.0026791|
| 2003-01-06 |    3|   0.0199587|   0.0130641|  -0.0016065|
| 2003-01-07 |    4|  -0.0133226|  -0.0074444|  -0.0113180|
| 2003-01-08 |    5|  -0.0270011|  -0.0188565|  -0.0169449|
| 2003-01-09 |    6|   0.0116346|   0.0169868|   0.0136876|

Ease computation by first standardizing the time series

``` r
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
                iterations = 100,
                P = 1, Q = 1,
                distribution = "Student_t",
                meanstructure = "constant")
#> 
#> CHECKING DATA AND PREPROCESSING FOR MODEL 'CCCMGARCH' NOW.
#> 
#> COMPILING MODEL 'CCCMGARCH' NOW.
#> 
#> STARTING SAMPLER FOR MODEL 'CCCMGARCH' NOW.
```

### Parameter Estimates

``` r
summary( fit1 )
#> Model: CCC-MGARCH
#> Basic Specification: H_t = D_t R D_t
#>  diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1
#> 
#> Distribution:  Student_t
#> ---
#> Iterations:  100
#> Chains:  4
#> Date:  Fri Aug 21 17:49:27 2020
#> Elapsed time (min):  0.37
#> 
#> GARCH(1,1)  estimates for conditional variance:
#> 
#>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> a_h_1,ty   0.10 0.07 0.08 0.01  0.27 319.87 0.99
#> a_h_1,ns   0.08 0.07 0.06 0.00  0.24 421.53 0.98
#> a_h_1,hn   0.11 0.08 0.08 0.01  0.29 460.21 0.98
#> b_h_1,ty   0.46 0.18 0.48 0.08  0.76 404.30 0.99
#> b_h_1,ns   0.39 0.18 0.37 0.10  0.79 343.41 1.00
#> b_h_1,hn   0.38 0.20 0.38 0.06  0.78 460.21 0.99
#> c_h_var_ty 0.29 0.12 0.27 0.11  0.55 386.13 1.00
#> c_h_var_ns 0.34 0.12 0.33 0.10  0.57 271.76 1.00
#> c_h_var_hn 0.46 0.19 0.43 0.16  0.87 460.21 0.98
#> 
#> 
#> Constant correlation (R) coefficients:
#> 
#>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> R_ns-ty 0.64 0.06 0.63 0.50  0.76 177.53 1.00
#> R_hn-ty 0.73 0.05 0.73 0.63  0.82 264.09 1.00
#> R_hn-ns 0.63 0.07 0.63 0.49  0.74 220.31 0.99
#> 
#> 
#> Intercept estimates on the location:
#> 
#>                     mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#> (Intercept)_toyota -0.09 0.08 -0.10 -0.23  0.07 185.50 1.00
#> (Intercept)_nissan -0.01 0.08 -0.01 -0.16  0.15 137.92 1.01
#> (Intercept)_honda  -0.03 0.09 -0.04 -0.20  0.15 166.08 1.00
#> 
#> 
#> Df constant student_t (nu):
#> 
#>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
#>  35.38  29.11  25.50   7.37 123.97 231.94   1.01 
#> 
#> 
#> Log density posterior estimate:
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#> -178.69    5.13 -178.04 -189.53 -170.11   55.13    1.09
```

### Forecasted Values

Forecast volatility 10 days ahead

``` r
fc <- forecast(fit1, ahead = 10 )
fc
#> ---
#> [Variance] Forecast for 10 ahead:
#> 
#> toyota :
#>       
#> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#>    101 0.54 0.11 0.53 0.36  0.76 259.69 1.00
#>    102 0.60 0.20 0.57 0.36  1.06 145.45 1.00
#>    103 0.62 0.25 0.58 0.35  1.09 190.82 1.00
#>    104 0.63 0.22 0.59 0.35  1.14 220.65 1.00
#>    105 0.62 0.17 0.60 0.38  1.05 202.52 1.00
#>    106 0.65 0.26 0.60 0.38  1.17 236.12 1.00
#>    107 0.66 0.25 0.62 0.38  1.33 220.28 1.01
#>    108 0.67 0.26 0.61 0.39  1.14 202.25 1.01
#>    109 0.67 0.29 0.60 0.38  1.30 185.98 1.00
#>    110 0.66 0.23 0.62 0.41  1.13 199.89 1.00
#> nissan :
#>       
#> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#>    101 0.60 0.11 0.60 0.44  0.87 221.34    1
#>    102 0.63 0.13 0.63 0.41  0.93 183.96    1
#>    103 0.64 0.16 0.62 0.40  0.93 208.40    1
#>    104 0.69 0.52 0.62 0.42  1.00 195.82    1
#>    105 0.66 0.32 0.62 0.41  0.95 192.86    1
#>    106 0.68 0.51 0.63 0.43  1.00 200.18    1
#>    107 0.65 0.21 0.63 0.42  1.01 188.26    1
#>    108 0.65 0.16 0.63 0.42  1.07 206.07    1
#>    109 0.67 0.20 0.64 0.42  1.14 208.63    1
#>    110 0.67 0.23 0.63 0.44  1.19 195.02    1
#> honda :
#>       
#> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#>    101 0.77 0.15 0.76 0.52  1.08 206.00 1.00
#>    102 0.84 0.21 0.79 0.55  1.35 184.63 0.99
#>    103 0.87 0.28 0.82 0.54  1.39 183.93 1.00
#>    104 0.90 0.43 0.84 0.53  1.68 200.65 1.00
#>    105 0.92 0.42 0.84 0.54  1.63 206.41 1.00
#>    106 0.92 0.34 0.84 0.56  1.80 222.99 1.00
#>    107 0.92 0.32 0.84 0.56  1.72 211.54 1.00
#>    108 0.94 0.36 0.86 0.57  1.68 217.35 1.00
#>    109 0.89 0.24 0.84 0.55  1.52 210.82 1.00
#>    110 0.90 0.25 0.84 0.55  1.43 211.51 1.00
```

``` r
plot(fc,askNewPage = FALSE, type = 'var' )
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-7-3.png" width="100%" />

Add two additional models, one with CCC(2,2) and a DCC(1,1)

``` r
# Fit CCC(1, 1) with constant on the mean structure.
fit2 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                parameterization = "CCC",
                iterations = 100,
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
                iterations = 100,
                P = 1, Q = 1,
                distribution = "Student_t",
                meanstructure = "arma")
#> 
#> CHECKING DATA AND PREPROCESSING FOR MODEL 'DCCMGARCH' NOW.
#> 
#> COMPILING MODEL 'DCCMGARCH' NOW.
#> 
#> STARTING SAMPLER FOR MODEL 'DCCMGARCH' NOW.
```

The DCC(1,1) model also incorportes an ARMA(1,1) meanstructure. The output will have the according information:

``` r
summary( fit3 )
#> Model: DCC-MGARCH
#> Basic Specification: H_t = D_t R D_t
#>  diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]
#>  R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])
#> 
#> Distribution:  Student_t
#> ---
#> Iterations:  100
#> Chains:  4
#> Date:  Fri Aug 21 17:51:45 2020
#> Elapsed time (min):  1.52
#> 
#> GARCH(1,1)  estimates for conditional variance on D:
#> 
#>            mean   sd  mdn 2.5% 97.5% n_eff Rhat
#> a_h_1,ty   0.16 0.14 0.13 0.01  0.56 34.99 1.06
#> a_h_1,ns   0.09 0.08 0.06 0.00  0.30 20.38 1.23
#> a_h_1,hn   0.10 0.08 0.08 0.00  0.28 25.40 1.18
#> b_h_1,ty   0.45 0.15 0.46 0.18  0.71 69.29 1.02
#> b_h_1,ns   0.36 0.23 0.36 0.01  0.80  4.22 1.43
#> b_h_1,hn   0.54 0.18 0.55 0.17  0.86 15.69 1.20
#> c_h_var_ty 0.27 0.11 0.26 0.12  0.49 50.26 1.10
#> c_h_var_ns 0.35 0.14 0.35 0.11  0.63  5.17 1.32
#> c_h_var_hn 0.32 0.13 0.30 0.12  0.61 38.98 1.11
#> 
#> 
#> GARCH(1,1) estimates for conditional variance on Q:
#> 
#>     mean   sd  mdn 2.5% 97.5% n_eff Rhat
#> a_q 0.18 0.10 0.17 0.02  0.39 19.72 1.18
#> b_q 0.26 0.19 0.23 0.03  0.69 16.44 1.23
#> 
#> 
#> Unconditional correlation 'S' in Q:
#> 
#>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> S_ns-ty 0.62 0.08 0.62 0.44  0.74 184.74 1.01
#> S_hn-ty 0.73 0.06 0.74 0.60  0.83 180.47 1.00
#> S_hn-ns 0.63 0.08 0.64 0.40  0.75 108.73 1.04
#> 
#> 
#> ARMA(1,1) estimates on the location:
#> 
#>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#> (Intercept)_toyota  -0.08 0.09 -0.08 -0.26  0.13 325.67 0.98
#> (Intercept)_nissan   0.02 0.09  0.00 -0.13  0.21 184.71 1.00
#> (Intercept)_honda   -0.02 0.12 -0.02 -0.27  0.23 260.86 1.00
#> Phi_toyota-toyota   -0.08 0.30 -0.12 -0.59  0.56  27.30 1.10
#> Phi_toyota-nissan    0.26 0.44  0.23 -0.55  0.95  12.62 1.35
#> Phi_toyota-honda     0.06 0.29  0.06 -0.43  0.59  26.89 1.10
#> Phi_nissan-toyota    0.15 0.36  0.17 -0.53  0.78  21.88 1.14
#> Phi_nissan-nissan   -0.05 0.37 -0.03 -0.80  0.70  31.99 1.11
#> Phi_nissan-honda     0.12 0.32  0.17 -0.50  0.68  53.06 1.05
#> Phi_honda-toyota    -0.26 0.36 -0.32 -0.90  0.43  34.57 1.11
#> Phi_honda-nissan     0.16 0.38  0.21 -0.77  0.80  45.71 1.13
#> Phi_honda-honda     -0.21 0.29 -0.20 -0.74  0.31  28.63 1.12
#> Theta_toyota-toyota -0.02 0.33  0.01 -0.79  0.53  21.37 1.14
#> Theta_toyota-nissan -0.10 0.42 -0.08 -0.83  0.70  13.13 1.31
#> Theta_toyota-honda  -0.03 0.30 -0.06 -0.53  0.59  28.93 1.10
#> Theta_nissan-toyota -0.18 0.36 -0.21 -0.84  0.52  21.48 1.17
#> Theta_nissan-nissan  0.08 0.33  0.09 -0.59  0.73  36.92 1.08
#> Theta_nissan-honda  -0.15 0.32 -0.20 -0.70  0.49  66.36 1.03
#> Theta_honda-toyota  -0.04 0.38  0.01 -0.74  0.73  33.21 1.12
#> Theta_honda-nissan  -0.03 0.41 -0.09 -0.71  0.91  37.57 1.17
#> Theta_honda-honda    0.36 0.32  0.38 -0.31  0.93  30.50 1.11
#> 
#> 
#> Df constant student_t (nu):
#> 
#>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
#>  45.50  28.61  39.17   7.99 107.14  89.08   1.02 
#> 
#> 
#> Log density posterior estimate:
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#> -176.38    6.01 -175.89 -189.77 -165.45   17.74    1.28
fc <- forecast(fit3, ahead =  10)
plot( fc,askNewPage = FALSE, type =  'mean' ) 
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-9-3.png" width="100%" />

### Ensemble Models

Obtain model weights with either the stacking or the pseudo BMA method. These methods are inherited from the `loo` package.

``` r
## use bmgarch_list function to collect bmgarch objects
modfits <- bmgarch_list(fit1, fit2, fit3)
```

### Compute Model Weights

Compute model weights with the stacking method and the the approximate leave-future-out cross validation (LFO CV). `L` defines the minimal length of the time series before we start engaging in cross-validation. Note that the standard is to use the approximate `backward` method to CV as itresults in fewest refits. Exact CV is also available with `exact` but not encouraged as it results in refitting all CV models.

``` r
mw <- model_weights(modfits, L = 50, method = 'stacking' )
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
#> Using threshold  0.6 , model was refit  4  times, at observations 90 80 71 63 
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
#> Using threshold  0.6 , model was refit  4  times, at observations 83 73 72 60 
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
#> Using threshold  0.6 , model was refit  13  times, at observations 92 87 84 81 78 73 72 64 63 59 58 57 53

## Return model weights:
mw
#> Method: stacking
#> ------
#>        weight
#> model1 0.616 
#> model2 0.384 
#> model3 0.000
```

### Weighted Forecasting

Use model weights to obtain weighted forecasts. Here we will forecast 10 days ahead.

``` r
w_fc <- forecast(modfits, ahead = 10, weights = mw )

w_fc
#> ---
#> LFO-weighted forecasts across  3 models.
#> ---
#> [Mean] Forecast for 10 ahead:
#> 
#> toyota :
#>       
#> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
#>    101 -0.08 0.55 -0.09 -1.33  1.08    NA   NA
#>    102 -0.11 0.56 -0.06 -1.56  0.85    NA   NA
#>    103 -0.09 0.62 -0.06 -1.52  1.04    NA   NA
#>    104 -0.13 0.59 -0.11 -1.27  0.97    NA   NA
#>    105 -0.04 0.67 -0.02 -1.61  1.21    NA   NA
#>    106 -0.05 0.64 -0.01 -1.18  1.19    NA   NA
#>    107 -0.09 0.64 -0.10 -1.22  1.20    NA   NA
#>    108 -0.03 0.67 -0.03 -1.33  1.23    NA   NA
#>    109 -0.13 0.66 -0.11 -1.31  1.23    NA   NA
#>    110 -0.01 0.62 -0.08 -1.20  1.23    NA   NA
#> nissan :
#>       
#> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
#>    101 -0.02 0.59 -0.01 -1.22  0.95    NA   NA
#>    102 -0.06 0.60 -0.02 -1.35  1.05    NA   NA
#>    103 -0.02 0.64  0.05 -1.53  1.05    NA   NA
#>    104 -0.07 0.60 -0.09 -1.10  1.24    NA   NA
#>    105  0.05 0.68  0.01 -1.25  1.45    NA   NA
#>    106 -0.01 0.62  0.00 -1.21  1.24    NA   NA
#>    107  0.01 0.57  0.02 -1.12  1.25    NA   NA
#>    108  0.06 0.66  0.06 -1.30  1.20    NA   NA
#>    109 -0.03 0.70  0.00 -1.37  1.24    NA   NA
#>    110  0.09 0.56  0.11 -1.12  1.15    NA   NA
#> honda :
#>       
#> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
#>    101 -0.04 0.64 -0.06 -1.25  1.04    NA   NA
#>    102 -0.05 0.70  0.04 -1.54  1.06    NA   NA
#>    103 -0.05 0.70 -0.01 -1.28  1.19    NA   NA
#>    104 -0.08 0.74 -0.06 -1.37  1.31    NA   NA
#>    105 -0.05 0.73  0.02 -1.49  1.16    NA   NA
#>    106  0.03 0.70 -0.06 -1.24  1.45    NA   NA
#>    107 -0.03 0.76 -0.09 -1.46  1.43    NA   NA
#>    108  0.04 0.69  0.04 -1.31  1.32    NA   NA
#>    109 -0.07 0.76 -0.04 -1.55  1.46    NA   NA
#>    110  0.07 0.67  0.04 -1.24  1.45    NA   NA
#> ---
#> [Variance] Forecast for 10 ahead:
#> 
#> toyota :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.53 0.08 0.53 0.38  0.68    NA   NA
#>    102 0.57 0.14 0.56 0.38  0.87    NA   NA
#>    103 0.60 0.17 0.57 0.38  0.94    NA   NA
#>    104 0.62 0.16 0.59 0.40  1.04    NA   NA
#>    105 0.62 0.14 0.60 0.39  0.93    NA   NA
#>    106 0.64 0.20 0.61 0.40  1.25    NA   NA
#>    107 0.66 0.24 0.62 0.42  1.15    NA   NA
#>    108 0.67 0.28 0.62 0.43  1.14    NA   NA
#>    109 0.68 0.26 0.62 0.43  1.24    NA   NA
#>    110 0.67 0.20 0.63 0.44  1.14    NA   NA
#> nissan :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.61 0.08 0.61 0.47  0.77    NA   NA
#>    102 0.63 0.10 0.62 0.46  0.86    NA   NA
#>    103 0.64 0.12 0.62 0.46  0.86    NA   NA
#>    104 0.68 0.33 0.62 0.49  1.01    NA   NA
#>    105 0.66 0.22 0.63 0.47  0.99    NA   NA
#>    106 0.68 0.32 0.64 0.50  1.02    NA   NA
#>    107 0.67 0.19 0.64 0.47  0.95    NA   NA
#>    108 0.66 0.13 0.65 0.45  0.97    NA   NA
#>    109 0.67 0.17 0.65 0.46  1.16    NA   NA
#>    110 0.67 0.17 0.64 0.46  1.00    NA   NA
#> honda :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.77 0.10 0.77 0.59  0.97    NA   NA
#>    102 0.81 0.14 0.80 0.60  1.15    NA   NA
#>    103 0.86 0.19 0.84 0.59  1.18    NA   NA
#>    104 0.88 0.28 0.85 0.59  1.41    NA   NA
#>    105 0.91 0.29 0.86 0.60  1.50    NA   NA
#>    106 0.91 0.29 0.85 0.61  1.65    NA   NA
#>    107 0.93 0.31 0.85 0.60  1.60    NA   NA
#>    108 0.92 0.27 0.86 0.63  1.48    NA   NA
#>    109 0.89 0.19 0.87 0.62  1.37    NA   NA
#>    110 0.91 0.21 0.88 0.62  1.32    NA   NA
#> [Correlation] Forecast for 10 ahead:
#> 
#> nissan_toyota :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    102 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    103 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    104 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    105 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    106 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    107 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    108 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    109 0.64 0.05 0.64 0.54  0.73    NA   NA
#>    110 0.64 0.05 0.64 0.54  0.73    NA   NA
#> honda_toyota :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    102 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    103 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    104 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    105 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    106 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    107 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    108 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    109 0.73 0.04 0.73 0.66   0.8    NA   NA
#>    110 0.73 0.04 0.73 0.66   0.8    NA   NA
#> honda_nissan :
#>       
#> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
#>    101 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    102 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    103 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    104 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    105 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    106 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    107 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    108 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    109 0.63 0.05 0.63 0.54  0.71    NA   NA
#>    110 0.63 0.05 0.63 0.54  0.71    NA   NA
plot(w_fc, askNewPage = FALSE, type =  'var' )
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-3.png" width="100%" />
