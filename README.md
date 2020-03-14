
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmgarch

`bmgarch` estimates Bayesian multivariate generalized autoregressive
conditional heteroskedasticity (MGARCH) models. Currently, bmgarch
supports ARMA(1,1) and intercept-only (Constant) mean structures, and a
variety of MGARCH(P,Q) parameterizations. In increasing order of
complexity:

  - CCC(P, Q): Constant Conditional Correlation
  - DCC(P, Q): Dynamic Conditional Correlation
  - BEKK(P, Q): Baba, Engle, Kraft, and Kroner
  - pdBEKK(P, Q): BEKK(P, Q) with positive diagonal constraints

## Installation

`bmgarch` is not yet available on CRAN.

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("ph-rast/bmgarch")
```

## Example

In this example, we use the pdBEKK(1,1) model for the variances, and an
intercept-only model for the means.

``` r
library(bmgarch)
#> Loading required package: Rcpp
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo

data(panas)
head(panas)
```

<div class="kable-table">

|     Pos |     Neg |
| ------: | ------: |
| \-2.193 | \-2.419 |
|   1.567 | \-0.360 |
| \-0.124 | \-1.202 |
|   0.020 | \-1.311 |
| \-0.150 |   2.004 |
|   3.877 |   1.008 |

</div>

``` r

# Fit pdBEKK(1, 1) with ARMA(1,1) on the mean structure.
fit <- bmgarch(panas,
               parameterization = "pdBEKK",
               iterations = 1000,
               P = 1, Q = 1,
               distribution = "Student_t",
               meanstructure = "constant")
#> 
#> CHECKING DATA AND PREPROCESSING FOR MODEL 'pdBEKKMGARCH' NOW.
#> 
#> COMPILING MODEL 'pdBEKKMGARCH' NOW.
#> 
#> STARTING SAMPLER FOR MODEL 'pdBEKKMGARCH' NOW.
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
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
#> Date:  Fri Mar 13 17:56:18 2020
#> Elapsed time (min):  5.77
#> 
#> ---
#> Constant correlation, R (diag[C]*R*diag[C]):
#> 
#>         mean  sd mdn  2.5% 97.5%   n_eff Rhat
#> R_Ng-Ps 0.01 0.5   0 -0.91  0.91 1657.52    1
#> 
#> 
#> Constant variances (diag[C]):
#> 
#>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> var_Ps 0.67 0.86 0.30 0.01  3.21 545.35 1.01
#> var_Ng 1.27 0.41 1.31 0.23  1.96 777.96 1.01
#> 
#> 
#> MGARCH(1,1) estimates for A:
#> 
#>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
#> A_Ps-Ps 0.35 0.10 0.35  0.14  0.55 1925.96    1
#> A_Ng-Ps 0.05 0.08 0.05 -0.10  0.20 1721.95    1
#> A_Ps-Ng 0.06 0.14 0.07 -0.24  0.33 1520.21    1
#> A_Ng-Ng 0.39 0.12 0.40  0.15  0.61 1233.61    1
#> 
#> 
#> MGARCH(1,1) estimates for B:
#> 
#>          mean   sd   mdn  2.5% 97.5%   n_eff Rhat
#> B_Ps-Ps  0.74 0.23  0.84  0.09  0.95  472.28 1.01
#> B_Ng-Ps -0.07 0.16 -0.08 -0.43  0.34  454.75 1.02
#> B_Ps-Ng  0.24 0.40  0.25 -0.72  1.06  645.33 1.00
#> B_Ng-Ng  0.33 0.19  0.31  0.02  0.72 1330.75 1.00
#> 
#> 
#> Intercept estimates on the location:
#> 
#>                 mean   sd  mdn  2.5% 97.5%   n_eff Rhat
#> (Intercept)_Pos 0.01 0.13 0.01 -0.25  0.26 2498.10    1
#> (Intercept)_Neg 0.11 0.10 0.11 -0.08  0.30 2419.84    1
#> 
#> 
#> Df constant student_t (nu):
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#>   51.69   26.91   46.58   14.45  111.95 2232.11    1.00 
#> 
#> 
#> Log density posterior estimate:
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#> -802.98    4.62 -802.63 -813.23 -794.99  357.17    1.02
```

### Forecasted values

``` r
fit.fc <- forecast(fit, ahead = 5)

fit.fc
#> ---
#> [Variance] Forecast for 5 ahead:
#> 
#> Pos :
#>       
#> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
#>    201 7.62 2.38 7.38 4.05 12.78  921.44    1
#>    202 7.07 3.22 6.54 3.44 14.01 1310.85    1
#>    203 6.75 3.57 6.00 3.18 15.45 1607.94    1
#>    204 6.57 4.17 5.44 3.09 16.45 1947.53    1
#>    205 6.39 4.69 5.12 2.98 17.09 1907.63    1
#> Neg :
#>       
#> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
#>    201 1.90 0.35 1.86 1.34  2.72 1437.36    1
#>    202 2.17 0.73 2.00 1.38  4.04 1753.59    1
#>    203 2.25 1.02 2.01 1.40  4.49 1924.90    1
#>    204 2.28 1.37 2.03 1.41  4.61 2002.03    1
#>    205 2.26 1.21 2.01 1.38  4.61 1966.20    1
#> [Correlation] Forecast for 5 ahead:
#> 
#> Neg_Pos :
#>       
#> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
#>    201 -0.06 0.16 -0.06 -0.38  0.24 1126.10    1
#>    202 -0.06 0.20 -0.07 -0.43  0.34 1876.16    1
#>    203 -0.05 0.21 -0.06 -0.45  0.42 1676.73    1
#>    204 -0.04 0.20 -0.05 -0.44  0.40 1439.03    1
#>    205 -0.04 0.20 -0.04 -0.44  0.39 2027.60    1
```

``` r
plot(fit.fc, askNewPage = FALSE, type = "var")
```

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />
