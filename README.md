
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
#> Warning: There were 2 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
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
#> Date:  Fri Mar 13 19:11:42 2020
#> Elapsed time (min):  6.32
#> 
#> ---
#> Constant correlation, R (diag[C]*R*diag[C]):
#> 
#>         mean   sd   mdn  2.5% 97.5%   n_eff Rhat
#> R_Ng-Ps    0 0.49 -0.01 -0.91  0.89 1425.38    1
#> 
#> 
#> Constant variances (diag[C]):
#> 
#>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
#> var_Ps 0.73 0.92 0.35 0.01  3.40 379.60 1.01
#> var_Ng 1.26 0.42 1.32 0.24  1.94 518.92 1.00
#> 
#> 
#> MGARCH(1,1) estimates for A:
#> 
#>         mean   sd  mdn  2.5% 97.5%   n_eff Rhat
#> A_Ps-Ps 0.34 0.11 0.35  0.11  0.55  912.66    1
#> A_Ng-Ps 0.06 0.08 0.05 -0.09  0.21 1642.67    1
#> A_Ps-Ng 0.06 0.14 0.06 -0.23  0.33 1488.82    1
#> A_Ng-Ng 0.39 0.12 0.39  0.13  0.60  825.11    1
#> 
#> 
#> MGARCH(1,1) estimates for B:
#> 
#>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
#> B_Ps-Ps  0.73 0.24  0.83  0.06  0.94 293.45 1.01
#> B_Ng-Ps -0.08 0.16 -0.08 -0.47  0.29 312.48 1.01
#> B_Ps-Ng  0.27 0.41  0.27 -0.73  1.13 520.00 1.00
#> B_Ng-Ng  0.33 0.20  0.32  0.01  0.74 933.35 1.00
#> 
#> 
#> Intercept estimates on the location:
#> 
#>                 mean   sd  mdn  2.5% 97.5%   n_eff Rhat
#> (Intercept)_Pos 0.01 0.14 0.00 -0.27   0.3 2349.14    1
#> (Intercept)_Neg 0.11 0.10 0.11 -0.09   0.3 2504.08    1
#> 
#> 
#> Df constant student_t (nu):
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#>   53.85   28.13   48.41   15.39  119.21 2180.27    1.00 
#> 
#> 
#> Log density posterior estimate:
#> 
#>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
#> -803.37    4.81 -802.97 -814.24 -795.31  330.79    1.02
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
#>    201 7.48 2.41 7.26 3.98 12.70  491.19    1
#>    202 7.01 3.12 6.51 3.41 14.59  916.84    1
#>    203 6.82 4.17 5.85 3.22 17.19 1624.79    1
#>    204 6.50 3.98 5.39 3.07 18.03 1543.71    1
#>    205 6.37 4.39 5.03 2.88 17.63 1821.36    1
#> Neg :
#>       
#> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
#>    201 1.91 0.35 1.86 1.35  2.78 1142.94    1
#>    202 2.17 0.81 1.99 1.39  4.01 1962.46    1
#>    203 2.25 1.28 2.01 1.37  4.34 1915.92    1
#>    204 2.27 1.04 2.03 1.39  4.70 1918.36    1
#>    205 2.27 1.08 2.00 1.39  5.15 2101.68    1
#> [Correlation] Forecast for 5 ahead:
#> 
#> Neg_Pos :
#>       
#> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
#>    201 -0.05 0.15 -0.05 -0.36  0.24  958.56    1
#>    202 -0.06 0.20 -0.06 -0.45  0.37 1352.61    1
#>    203 -0.05 0.20 -0.05 -0.45  0.37 1544.84    1
#>    204 -0.04 0.21 -0.05 -0.43  0.42 1537.14    1
#>    205 -0.02 0.20 -0.02 -0.43  0.41 2046.01    1
```

``` r
plot(fit.fc, askNewPage = FALSE, type = "var")
```

<img src="man/figures/README-forecastPlot-1.png" width="100%" /><img src="man/figures/README-forecastPlot-2.png" width="100%" />

``` r

plot(fit.fc, askNewPage = FALSE, type = "cor")
```

<img src="man/figures/README-forecastPlot-3.png" width="100%" />
