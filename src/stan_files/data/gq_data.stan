int<lower=2> T;
int<lower=2> nt;    // number of time series
int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
int<lower=1> P; // AR component in MGARCH(P,Q), matrix B  
vector[nt] rts[T];  // multivariate time-series
vector[nt] xC[T];  // time-varying predictor for conditional H
int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
int<lower=0, upper=1> meanstructure; // Select model for location
int<lower=1> ahead; // forecasted periods
vector[nt] xC_p[ahead];  // time-varying predictor for conditional H
vector[nt] future_rts[ahead]; // future observations to obtain log_lik
int<lower = 0, upper = 1> compute_log_lik;
