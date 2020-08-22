<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with render("README.Rmd", output_format = "md_document") -->
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
    #> Date:  Sat Aug 22 13:00:18 2020
    #> Elapsed time (min):  15.45
    #> 
    #> ---
    #> Constant correlation, R (diag[C]*R*diag[C]):
    #> 
    #>         mean   sd  mdn  2.5% 97.5% n_eff Rhat
    #> R_Ng-Ps  0.2 0.56 0.24 -0.89  0.93  4.96  1.3
    #> 
    #> 
    #> Constant variances (diag[C]):
    #> 
    #>        mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> var_Ps 0.50 0.78 0.16 0.01  3.01  21.55 1.10
    #> var_Ng 1.19 0.37 1.20 0.34  1.88 486.58 1.01
    #> 
    #> 
    #> MGARCH(1,1) estimates for A:
    #> 
    #>         mean   sd  mdn  2.5% 97.5% n_eff Rhat
    #> A_Ps-Ps 0.36 0.11 0.36  0.13  0.55 18.83 1.14
    #> A_Ng-Ps 0.07 0.07 0.07 -0.08  0.22 38.57 1.08
    #> A_Ps-Ng 0.01 0.14 0.01 -0.25  0.30 25.42 1.11
    #> A_Ng-Ng 0.42 0.13 0.42  0.13  0.64 14.48 1.13
    #> 
    #> 
    #> MGARCH(1,1) estimates for B:
    #> 
    #>          mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> B_Ps-Ps  0.76 0.22  0.84  0.12  0.94  58.75 1.07
    #> B_Ng-Ps -0.10 0.15 -0.10 -0.46  0.23 348.80 1.01
    #> B_Ps-Ng  0.32 0.37  0.34 -0.61  1.13 290.03 1.02
    #> B_Ng-Ng  0.37 0.18  0.38  0.03  0.75 267.70 1.02
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_Pos -0.08 0.17 -0.07 -0.38  0.26   7.47 1.20
    #> (Intercept)_Neg  0.05 0.12  0.03 -0.15  0.34  20.94 1.09
    #> Phi_Pos-Pos     -0.15 0.37 -0.17 -0.78  0.59   7.56 1.20
    #> Phi_Pos-Neg     -0.28 0.51 -0.36 -0.92  0.75   5.69 1.27
    #> Phi_Neg-Pos     -0.16 0.32 -0.17 -0.78  0.53 204.91 1.03
    #> Phi_Neg-Neg      0.16 0.43  0.24 -0.75  0.80   7.91 1.19
    #> Theta_Pos-Pos    0.06 0.39  0.11 -0.69  0.73   8.55 1.18
    #> Theta_Pos-Neg    0.20 0.53  0.26 -0.87  0.86   5.13 1.31
    #> Theta_Neg-Pos    0.19 0.32  0.20 -0.51  0.80 192.64 1.04
    #> Theta_Neg-Neg   -0.19 0.45 -0.26 -0.84  0.76   7.38 1.21
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  53.82  24.84  51.43  14.54 110.60 344.10   1.02 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -802.89    4.58 -802.44 -813.85 -794.95  215.88    1.04

### Forecasted values

    fit.fc <- forecast(fit, ahead = 5)

    fit.fc
    #> ---
    #> [Mean] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%   n_eff Rhat
    #>    201 -0.44 3.01 -0.38 -6.44  5.53 1981.31    1
    #>    202 -0.32 3.04 -0.31 -6.30  5.56 1821.94    1
    #>    203 -0.24 2.96 -0.27 -5.96  5.42 1838.34    1
    #>    204 -0.18 2.74 -0.09 -5.77  5.17 1968.75    1
    #>    205 -0.15 2.75 -0.19 -5.73  5.52 1897.94    1
    #> Neg :
    #>       
    #> period mean   sd  mdn  2.5% 97.5%   n_eff Rhat
    #>    201 0.48 1.51 0.44 -2.49  3.36 1267.81 1.01
    #>    202 0.32 1.61 0.30 -2.92  3.58 1110.87 1.01
    #>    203 0.28 1.60 0.25 -2.84  3.45 1793.79 1.00
    #>    204 0.19 1.61 0.18 -2.98  3.46 1631.36 1.01
    #>    205 0.20 1.62 0.16 -3.03  3.31 2019.79 1.00
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> Pos :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #>    201 8.96 3.15 8.84 4.15 15.13  16.40 1.19
    #>    202 8.16 4.43 7.47 3.46 18.41  61.35 1.07
    #>    203 7.69 5.45 6.48 3.28 19.55 184.91 1.04
    #>    204 7.38 6.37 5.82 3.06 22.00 444.59 1.02
    #>    205 7.07 6.02 5.38 2.82 22.42 504.82 1.01
    #> Neg :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%   n_eff Rhat
    #>    201 1.95 0.41 1.89 1.35  2.91  198.14 1.00
    #>    202 2.27 0.87 2.05 1.40  4.32  638.74 1.00
    #>    203 2.38 1.12 2.10 1.41  5.55  893.52 1.00
    #>    204 2.40 1.46 2.08 1.40  5.12 1177.66 1.00
    #>    205 2.40 1.28 2.03 1.38  5.83  721.93 1.01
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> Neg_Pos :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #>    201 -0.11 0.17 -0.11 -0.43  0.21  39.73 1.07
    #>    202 -0.10 0.22 -0.11 -0.51  0.34 170.46 1.03
    #>    203 -0.08 0.23 -0.07 -0.51  0.41 454.97 1.02
    #>    204 -0.05 0.23 -0.05 -0.48  0.44 629.26 1.01
    #>    205 -0.04 0.22 -0.03 -0.49  0.42 878.96 1.00

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
                    iterations = 100,
                    P = 1, Q = 1,
                    distribution = "Student_t",
                    meanstructure = "constant")

### Parameter Estimates

    summary( fit1 )
    #> Model: CCC-MGARCH
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1
    #> 
    #> Distribution:  Student_t
    #> ---
    #> Iterations:  100
    #> Chains:  4
    #> Date:  Sat Aug 22 13:00:52 2020
    #> Elapsed time (min):  0.31
    #> 
    #> GARCH(1,1)  estimates for conditional variance:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_h_1,ty   0.10 0.09 0.07 0.00  0.32 216.95 0.99
    #> a_h_1,ns   0.07 0.06 0.06 0.00  0.22 318.50 1.00
    #> a_h_1,hn   0.11 0.08 0.10 0.01  0.32 266.14 1.00
    #> b_h_1,ty   0.46 0.19 0.47 0.10  0.78 345.05 1.00
    #> b_h_1,ns   0.38 0.21 0.37 0.06  0.80 460.21 0.99
    #> b_h_1,hn   0.40 0.17 0.39 0.11  0.73 338.31 0.99
    #> c_h_var_ty 0.29 0.13 0.26 0.09  0.58 258.48 1.00
    #> c_h_var_ns 0.36 0.13 0.35 0.11  0.59 460.21 0.98
    #> c_h_var_hn 0.44 0.16 0.42 0.13  0.78 324.13 0.99
    #> 
    #> 
    #> Constant correlation (R) coefficients:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> R_ns-ty 0.66 0.06 0.66 0.52  0.77 189.00 0.99
    #> R_hn-ty 0.73 0.05 0.74 0.63  0.82 150.61 1.00
    #> R_hn-ns 0.64 0.06 0.65 0.53  0.76 185.10 0.99
    #> 
    #> 
    #> Intercept estimates on the location:
    #> 
    #>                     mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota -0.09 0.07 -0.08 -0.23  0.05 131.05 1.00
    #> (Intercept)_nissan  0.00 0.08  0.00 -0.15  0.13 154.58 1.01
    #> (Intercept)_honda  -0.03 0.08 -0.03 -0.20  0.12 167.42 1.00
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  33.50  23.17  26.72   7.92  88.88 317.41   0.99 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -178.20    5.72 -177.83 -190.31 -167.96   65.99    1.02

### Forecasted Values

Forecast volatility 10 days ahead

    fc <- forecast(fit1, ahead = 10 )
    fc
    #> ---
    #> [Variance] Forecast for 10 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #>    101 0.55 0.11 0.55 0.36  0.77 179.87 1.01
    #>    102 0.59 0.16 0.58 0.36  0.86 176.19 1.00
    #>    103 0.63 0.19 0.60 0.41  1.08 220.52 1.00
    #>    104 0.65 0.26 0.60 0.41  1.27 200.75 1.00
    #>    105 0.65 0.23 0.60 0.42  1.12 195.19 1.00
    #>    106 0.71 0.69 0.62 0.40  1.37 218.57 1.00
    #>    107 0.77 1.37 0.63 0.41  1.39 209.22 1.00
    #>    108 0.75 1.17 0.62 0.40  1.64 209.84 1.00
    #>    109 0.77 1.19 0.62 0.41  1.43 214.09 1.00
    #>    110 0.72 0.50 0.63 0.40  1.47 225.24 1.01
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #>    101 0.62 0.11 0.60 0.44  0.86 167.80 1.00
    #>    102 0.65 0.14 0.63 0.45  0.95 166.31 1.00
    #>    103 0.66 0.16 0.63 0.45  1.07 146.21 1.00
    #>    104 0.67 0.22 0.63 0.44  1.09 203.13 1.00
    #>    105 0.67 0.16 0.64 0.45  1.01 218.14 1.00
    #>    106 0.66 0.14 0.64 0.45  0.99 212.17 1.00
    #>    107 0.67 0.21 0.63 0.45  1.06 191.47 1.00
    #>    108 0.68 0.22 0.65 0.44  1.06 211.28 1.01
    #>    109 0.70 0.38 0.64 0.44  1.28 217.78 1.01
    #>    110 0.69 0.24 0.65 0.45  1.06 177.89 1.01
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #>    101 0.77 0.13 0.77 0.54  1.01 155.89 1.01
    #>    102 0.84 0.17 0.83 0.57  1.32 156.61 1.00
    #>    103 0.86 0.21 0.84 0.56  1.39 188.84 1.00
    #>    104 0.90 0.31 0.84 0.59  1.71 198.56 1.00
    #>    105 0.92 0.36 0.86 0.56  1.60 238.64 1.00
    #>    106 0.91 0.30 0.85 0.56  1.57 199.05 1.00
    #>    107 0.92 0.32 0.85 0.56  1.82 186.78 1.00
    #>    108 0.91 0.28 0.87 0.54  1.74 146.75 1.00
    #>    109 0.99 0.64 0.86 0.57  2.02 180.71 1.00
    #>    110 0.95 0.44 0.87 0.55  1.87 165.26 1.00

    plot(fc,askNewPage = FALSE, type = 'var' )

<img src="man/figures/README-stockForecastPlot-1.png" width="100%" /><img src="man/figures/README-stockForecastPlot-2.png" width="100%" /><img src="man/figures/README-stockForecastPlot-3.png" width="100%" />

Add two additional models, one with CCC(2,2) and a DCC(1,1)

    # Fit CCC(1, 1) with constant on the mean structure.
    fit2 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                    parameterization = "CCC",
                    iterations = 100,
                    P = 2, Q = 2,
                    distribution = "Student_t",
                    meanstructure = "constant")

    fit3 <- bmgarch(stocks.z[1:100, c("toyota", "nissan", "honda")],
                    parameterization = "DCC",
                    iterations = 100,
                    P = 1, Q = 1,
                    distribution = "Student_t",
                    meanstructure = "arma")

The DCC(1,1) model also incorportes an ARMA(1,1) meanstructure. The
output will have the according information:

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
    #> Date:  Sat Aug 22 13:03:08 2020
    #> Elapsed time (min):  1.49
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_h_1,ty   0.16 0.12 0.13 0.02  0.44 130.05 1.03
    #> a_h_1,ns   0.09 0.08 0.07 0.00  0.29 193.76 1.01
    #> a_h_1,hn   0.12 0.10 0.10 0.01  0.36  25.22 1.07
    #> b_h_1,ty   0.45 0.17 0.47 0.10  0.71  42.61 1.05
    #> b_h_1,ns   0.42 0.19 0.42 0.09  0.77 148.07 1.00
    #> b_h_1,hn   0.48 0.20 0.47 0.14  0.84  13.78 1.13
    #> c_h_var_ty 0.27 0.12 0.26 0.09  0.53  93.92 1.02
    #> c_h_var_ns 0.32 0.12 0.32 0.12  0.52 148.48 0.99
    #> c_h_var_hn 0.36 0.16 0.36 0.09  0.71  20.99 1.11
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>     mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> a_q 0.21 0.10 0.19 0.02  0.39 125.63 1.03
    #> b_q 0.21 0.14 0.21 0.01  0.49  27.08 1.08
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%  n_eff Rhat
    #> S_ns-ty 0.60 0.09 0.61 0.41  0.73 182.87 1.01
    #> S_hn-ty 0.74 0.06 0.74 0.58  0.84 206.03 0.99
    #> S_hn-ns 0.63 0.08 0.64 0.47  0.77 221.25 0.99
    #> 
    #> 
    #> ARMA(1,1) estimates on the location:
    #> 
    #>                      mean   sd   mdn  2.5% 97.5%  n_eff Rhat
    #> (Intercept)_toyota  -0.08 0.09 -0.07 -0.30  0.10 168.54 1.01
    #> (Intercept)_nissan   0.01 0.09  0.01 -0.18  0.18 259.92 0.99
    #> (Intercept)_honda    0.00 0.11  0.00 -0.18  0.22 143.61 1.01
    #> Phi_toyota-toyota    0.01 0.34  0.01 -0.62  0.67 102.59 1.02
    #> Phi_toyota-nissan   -0.11 0.40 -0.16 -0.78  0.70  15.77 1.18
    #> Phi_toyota-honda     0.14 0.40  0.12 -0.63  0.87  31.56 1.13
    #> Phi_nissan-toyota    0.24 0.38  0.31 -0.65  0.83  49.31 1.04
    #> Phi_nissan-nissan   -0.16 0.35 -0.19 -0.78  0.55  35.95 1.05
    #> Phi_nissan-honda     0.03 0.41  0.06 -0.91  0.79  15.44 1.15
    #> Phi_honda-toyota    -0.23 0.40 -0.25 -0.97  0.51  38.87 1.06
    #> Phi_honda-nissan     0.15 0.43  0.19 -0.74  0.89   9.90 1.18
    #> Phi_honda-honda     -0.10 0.34 -0.07 -0.75  0.49  81.00 1.00
    #> Theta_toyota-toyota -0.09 0.39 -0.12 -0.81  0.67  57.36 1.06
    #> Theta_toyota-nissan  0.24 0.40  0.27 -0.61  0.93  21.72 1.15
    #> Theta_toyota-honda  -0.13 0.38 -0.10 -0.83  0.56  39.80 1.10
    #> Theta_nissan-toyota -0.24 0.39 -0.29 -0.85  0.58  41.97 1.06
    #> Theta_nissan-nissan  0.16 0.34  0.19 -0.49  0.81  39.72 1.04
    #> Theta_nissan-honda  -0.08 0.39 -0.11 -0.78  0.71  27.42 1.09
    #> Theta_honda-toyota  -0.02 0.44 -0.04 -0.83  0.86  40.08 1.06
    #> Theta_honda-nissan  -0.06 0.46 -0.11 -0.80  0.86   8.07 1.21
    #> Theta_honda-honda    0.20 0.39  0.17 -0.49  0.94  82.27 1.01
    #> 
    #> 
    #> Df constant student_t (nu):
    #> 
    #>   mean     sd    mdn   2.5%  97.5%  n_eff   Rhat 
    #>  42.87  25.56  36.30  10.85 103.34 161.57   1.03 
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>    mean      sd     mdn    2.5%   97.5%   n_eff    Rhat 
    #> -176.60    5.35 -176.23 -187.29 -166.62   73.10    1.01
    fc <- forecast(fit3, ahead =  10)

    plot( fc,askNewPage = FALSE, type =  'mean' ) 

<img src="man/figures/README-fit3ForecastPlot-1.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-2.png" width="100%" /><img src="man/figures/README-fit3ForecastPlot-3.png" width="100%" />

### Ensemble Models

Obtain model weights with either the stacking or the pseudo BMA method.
These methods are inherited from the `loo` package.

    ## use bmgarch_list function to collect bmgarch objects
    modfits <- bmgarch_list(fit1, fit2, fit3)

### Compute Model Weights

Compute model weights with the stacking method and the the approximate
leave-future-out cross validation (LFO CV). `L` defines the minimal
length of the time series before we start engaging in cross-validation.
Note that the standard is to use the approximate `backward` method to CV
as itresults in fewest refits. Exact CV is also available with `exact`
but not encouraged as it results in refitting all CV models.

    mw <- model_weights(modfits, L = 50, method = 'stacking' )
    #> Using threshold  0.6 , model was refit  11  times, at observations 83 81 73 72 71 66 63 60 57 52 51

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
    #>    101 -0.17 0.80 -0.18 -1.96  1.37    NA   NA
    #>    102 -0.07 0.77 -0.12 -1.44  1.45    NA   NA
    #>    103 -0.11 0.78 -0.04 -1.77  1.41    NA   NA
    #>    104 -0.17 0.85 -0.17 -1.97  1.28    NA   NA
    #>    105 -0.14 0.84 -0.14 -1.97  1.33    NA   NA
    #> nissan :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101 -0.06 0.85 -0.08 -1.63  1.70    NA   NA
    #>    102  0.00 0.85 -0.01 -1.42  1.57    NA   NA
    #>    103  0.00 0.81  0.03 -1.71  1.46    NA   NA
    #>    104 -0.09 0.86  0.08 -1.99  1.53    NA   NA
    #>    105 -0.03 0.86 -0.02 -1.88  1.54    NA   NA
    #> honda :
    #>       
    #> period  mean   sd   mdn  2.5% 97.5% n_eff Rhat
    #>    101  0.00 0.94  0.03 -1.95  1.90    NA   NA
    #>    102 -0.08 0.90 -0.07 -1.85  1.81    NA   NA
    #>    103 -0.05 0.95 -0.02 -2.04  1.76    NA   NA
    #>    104 -0.20 1.00 -0.21 -2.47  1.65    NA   NA
    #>    105 -0.07 0.94 -0.05 -1.96  1.81    NA   NA
    #> ---
    #> [Variance] Forecast for 5 ahead:
    #> 
    #> toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.52 0.11 0.51 0.31  0.74    NA   NA
    #>    102 0.54 0.12 0.53 0.30  0.80    NA   NA
    #>    103 0.59 0.20 0.57 0.32  1.01    NA   NA
    #>    104 0.59 0.17 0.58 0.31  0.89    NA   NA
    #>    105 0.63 0.25 0.60 0.31  1.06    NA   NA
    #> nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.12 0.62 0.45  0.91    NA   NA
    #>    102 0.65 0.13 0.62 0.45  0.98    NA   NA
    #>    103 0.66 0.16 0.63 0.45  1.05    NA   NA
    #>    104 0.67 0.19 0.63 0.44  1.06    NA   NA
    #>    105 0.68 0.20 0.63 0.45  1.23    NA   NA
    #> honda :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.78 0.14 0.76 0.55  1.07    NA   NA
    #>    102 0.78 0.19 0.75 0.53  1.21    NA   NA
    #>    103 0.86 0.28 0.80 0.54  1.64    NA   NA
    #>    104 0.86 0.25 0.81 0.54  1.53    NA   NA
    #>    105 0.89 0.29 0.83 0.53  1.63    NA   NA
    #> [Correlation] Forecast for 5 ahead:
    #> 
    #> nissan_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.65 0.06 0.65 0.53  0.75    NA   NA
    #>    102 0.65 0.06 0.65 0.53  0.75    NA   NA
    #>    103 0.65 0.06 0.65 0.53  0.75    NA   NA
    #>    104 0.65 0.06 0.65 0.53  0.75    NA   NA
    #>    105 0.65 0.06 0.65 0.53  0.75    NA   NA
    #> honda_toyota :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.73 0.05 0.74 0.61  0.83    NA   NA
    #>    102 0.73 0.05 0.74 0.61  0.83    NA   NA
    #>    103 0.73 0.05 0.74 0.61  0.83    NA   NA
    #>    104 0.73 0.05 0.74 0.61  0.83    NA   NA
    #>    105 0.73 0.05 0.74 0.61  0.83    NA   NA
    #> honda_nissan :
    #>       
    #> period mean   sd  mdn 2.5% 97.5% n_eff Rhat
    #>    101 0.64 0.06 0.65 0.52  0.76    NA   NA
    #>    102 0.64 0.06 0.65 0.52  0.76    NA   NA
    #>    103 0.64 0.06 0.65 0.52  0.76    NA   NA
    #>    104 0.64 0.06 0.65 0.52  0.76    NA   NA
    #>    105 0.64 0.06 0.65 0.52  0.76    NA   NA

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
