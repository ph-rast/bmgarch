## Project title:
##    Created at:
##        Author: Philippe Rast
##          Data:
##       Summary:
## ---------------------------------------------------------------------- ##
rm(list=ls())
############################
## Read person level data ##
############################
ddip.names <-  scan(file = '~/UZH/Projekte/Archiv/LSM_Ferrer/WORK/-DATASETS/SOURCE/ddip_visits_Names.dat', nlines = 4, what = "character")
ddip.l23 <-  read.table(file = '~/UZH/Projekte/Archiv/LSM_Ferrer/WORK/-DATASETS/SOURCE/ddip_visits.dat', na.strings = '.',
                        col.names = ddip.names)

head(ddip.l23)
ddip.l23$id_individual


######################
## Data Wrangling   
##  Load ddip.lag
##  merge files
##  Select couples with at least 85 days
######################
load(file = "~/UZH/Projekte/M_GARCH/WORK/-DATASETS/DERIVED/ddip.rda")

## time series don't match up - some have provided  more datapoints than their partner
dyads <- unique(ddip.lag$ID_DYAD)
dyads

str(ddip.lag)

i=1
id <- unique(ddip.lag[ddip.lag$ID_DYAD==dyads[i], "ID_INDIV"])
out <- merge(x=ddip.lag[ddip.lag$ID_INDIV==id[1],], y=ddip.lag[ddip.lag$ID_INDIV==id[2],], by = 'DAY')
names(out)
dim(out)

names(ddip.lag)
header <- unlist(strsplit(names(out[2:30]), '.x'))
header
colnames(out[2:30]) <- header
names(out[31:59]) <- header
colnames(out) <- c('DAY', header, header)
out[1:10,]

out <- rbind(out[1:30],out[c(1, 31:59)]) 

for(i in 2:length(dyads)){
    id <- unique(ddip.lag[ddip.lag$ID_DYAD==dyads[i], "ID_INDIV"])
    out2 <- merge(x=ddip.lag[ddip.lag$ID_INDIV==id[1],], y=ddip.lag[ddip.lag$ID_INDIV==id[2],], by = 'DAY')
    colnames(out2) <- c('DAY', header, header)
    out2 <- rbind(out2[1:30],out2[c(1, 31:59)]) 
    out <- rbind(out, out2)
}
    
dim(out)
names(out)

ddip.lag <- out

names(ddip.lag)
names(ddip.l23)
ddip.full <- merge(x = ddip.lag, y = ddip.l23, by.x = "ID_INDIV", by.y = "id_individual", all.y = TRUE, incomparables = )
names(ddip.full)

ddip.full$ID_DYAD

ddip.full2 <- ddip.full[!is.na(ddip.full[,'ID_DYAD']),]

ddip.full2[1:150, 1:5]


## Selecton only those with at least 85 days
sel <-  unique(ddip.full2[ddip.full2$DAY>= 85, "ID_DYAD"])
sel

length(sel)

## Only select participants from sel
ddip.lag <- subset(ddip.full2, subset = ID_DYAD %in% sel)
names(ddip.lag)

## DAY does not mean that days are consecutive. Could  be: 1, 2, 5, 19, 85 
i_nobs_fem = NA
i_nobs_mal = NA
for(i in 1:length(sel) ){
    tst <- ddip.lag[ddip.lag$ID_DYAD == sel[i], 1:5]
    ind <- unique(tst$ID_INDIV)
    fem <- tst[tst$ID_INDIV==ind[1],]
    mal <- tst[tst$ID_INDIV==ind[2],]
    i_nobs_fem[i] = length(unique(sort(fem$DAY)))
    i_nobs_mal[i] = length(unique(sort(mal$DAY)))
}

cbind(i_nobs_fem, i_nobs_mal)

## drop unused dasets:
rm(list = c('ddip.full2', 'ddip.full', 'ddip.l23', 'tst'))



## some descriptive variables;
range(ddip.lag$V1age_yr)
mean(ddip.lag$V1age_yr)
sd(ddip.lag$V1age_yr)

range(ddip.lag$totyrinrel)
mean(ddip.lag$totyrinrel)
sd(ddip.lag$totyrinrel)

ddip.lag$pa_rel_std <- scale(ddip.lag$pa_rel)

#########################
## END Datawrangling   ##
#########################


##################################
## Load precompiled BEKK        ##
##################################


## Select one dyad
i <- 23
partner <- unique(ddip.lag[ddip.lag$ID_DYAD==sel[i], "ID_INDIV"])
partner
## Extract Time series per invididul in dyad
r1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1],'pa_rel_std']
r2 <- ddip.lag[ddip.lag$ID_INDIV==partner[2],'pa_rel_std']

pred <- ddip.lag[ddip.lag$ID_INDIV==partner[2],'na_rel']
## ## It seems some have almost no variability in NA (eg. i = 27)
## na1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1], 'na_rel']
## na2 <- ddip.lag[ddip.lag$ID_INDIV==partner[2], 'na_rel']

plot(r1, type = 'l', ylim = c(-3, 3), main = i)
lines(r2, col = 'red')

dev.off()

r <- cbind(r1, r2)
#r2 <- cbind(( (r - mean(c(r,pred))) / sd(r)), ( (pred - mean(c(r,pred))) / sd(r)))
sigma1 <- var(r)

r
#foreign::write.dta( as.data.frame(  r  ), file = '~/Downloads/dyad10.dta')

## Fit Model
fit <- bmgarch(data = r[, 1:2], xH = NULL,
               parameterization = "DCC", P = 1, Q = 1,
               iterations = 300,
               meanstructure = "constant")

summary(fit)

colnames(as.matrix(fit$model_fit))

plot(fit, type = "cvar")

aussi <-  forecast(fit , ahead =  2)
aussi

seed <-  NA


stanmodels$forecastDCC
 
str(aussi )


standata <-  list(T =  fit$TS_length,  nt =  fit$nt,
                 rts =  cbind( fit$RTS_full),
                 xH =  fit$xH,
                 Q =   fit$mgarchQ,
                 P =   fit$mgarchP,
                 ahead =  3, 
                 meanstructure =   fit$meanstructure,
                 distribution =   fit$num_dist)


as.matrix(fit$model_fit)[,'vd[2]']

frcst <- rstan::stan_model(file = "../../src/stan_files/forecastDCC.stan")

out <- rstan::gqs(frcst, draws = as.matrix(fit$model_fit), data =  standata)
out

dim(out)

str(out)

## return names and obtain position of rts_p's
pred_rtsp <- grep('rts_p', out@sim$fnames_oi)
pred_rtsp

out@sim$samples[[1]][pred_rtsp[1]]



rstan::summary(fit$model_fit, pars = c('vd'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')]

bekk_fit = bmgarch(data = r[,1:2], xH = NULL, parameterization = "BEKK", iterations = 500, P = 1, Q = 1,
                   meanstructure = 'constant')

bekk_fit = bmgarch(data = r[,1:2], xH = NULL, parameterization = "DCC", iterations = 500, P = 2, Q = 2,
                   meanstructure = 'arma')


summary(bekk_fit)

plot(bekk_fit, type = 'ccor')

pst <- out <- rstan::gqs(frcst, draws = as.matrix(bekk_fit$model_fit), data =  standat)

pst

## Not run:
library(rstan )

m <- stan_model(model_code ='
data {
  int<lower = 0, upper = 1> flag;
}
parameters {
  real y;
  vector[flag ? 1 : 0] phi;
}
transformed parameters {
  real mu;
  if ( flag == 1 ) {
    mu = 0.0 + phi[1];
  } else if ( flag == 0 ) {
    mu = 0.0;
  }  
}
model {
  phi ~ normal(10, 0.1);
  y ~ normal(mu, 1);
}')

f <- sampling(m, iter = 300, data =  list( flag =  1 ))

f

mc <-'
data {
  int<lower = 0, upper = 1> flag;
}
parameters {
  real y;
//  vector[flag ? 1 : 0] phi;
  real mu;
}
generated quantities {
//  real mu_p;
  real y_rep;
//  if ( flag == 1 ) {
//    mu_p = 0.0 + phi[1];
//  } else if ( flag == 0 ) {
//    mu_p = 0.0;
//  }
  y_rep = normal_rng(0, 1);
}'

mc <-'
data {
  int<lower = 0, upper = 1> flag;
}
parameters {
  real y;
//  vector[flag ? 1 : 0] phi;
  real mu;
}
generated quantities {
  real y_rep;
  y_rep = normal_rng(y, 1);
}'

m2 <- stan_model(model_code = mc)

f2 <- rstan::gqs(m2, draws = as.matrix(f), data =  list( flag =  0))
f2




###########################################
## Loop through participants
library(tcltk)

length(sel)
results_bekk <- array(NA, dim = c(length(sel),19))
resultsCI_bekk <- array(NA, dim = c(length(sel),38))

## loop through sel participants
total <- length(sel)
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

for(i in 1:length(sel)){
    
#for(i in redo){
    partner <- unique(ddip.lag[ddip.lag$ID_DYAD==sel[i], "ID_INDIV"])
    ## Extract Time series per invididul in dyad
    r1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1],'pa_rel_std']
    r2 <- ddip.lag[ddip.lag$ID_INDIV==partner[2],'pa_rel_std']
    r <-  cbind(r1, r2)
    fit1 = bmgarch(data = r, parameterization = "BEKK", iterations = 500, standardize_data = FALSE)

    pst <- rstan::summary(fit1$model_fit, probs = c(.05, .95), pars = c('Cnst', 'A', 'B', 'corC',  'phi0', 'phi', 'theta', 'lp__'))$summary[,c('mean', '5%', '95%', 'Rhat')]
  
    ## Overwrite parameter if Rhat < 1.1
    results_bekk[i, ] <- pst[,1]*ifelse(pst[,'Rhat']>1.1, NA, 1)
    colnames(results_bekk) <- rownames(pst)
    ## Record lower and upper CI
    resultsCI_bekk[i,] <-  c(t(pst[,2:3]))
    low <- paste0(rownames(pst), '5%')
    hi <- paste0(rownames(pst), '95%')
    colnames(resultsCI_bekk)  <-   c(t(cbind(low, hi)))    
    ## Progress Bar
    setTxtProgressBar(pb, i)
}
close(pb)

## focus on b2 for now?

results_bekk
resultsCI_bekk

redo <- which(is.na(rowSums(results_bekk)))
length(redo)
redo
total



save(results_bekk, file = 'results_dyad_bekk.RDat', compress = 'xz')
save(resultsCI_bekk, file = 'resultsCI_dyad_bekk.RDat', compress = 'xz')

load(file = 'results_dyad_bekk.RDat')
load(file = 'resultsCI_dyad_bekk.RDat')
    
## Create filter for results with posterior mass outside of zero
    nonzero <- ifelse(sign(resultsCI_bekk[, seq(1, ncol(resultsCI_bekk), by = 2)]) +
                      sign(resultsCI_bekk[, seq(2, ncol(resultsCI_bekk), by = 2)]) == 0, NA, 1)

## Return nonzero results:
round(results_bekk * nonzero, 2)

sel[22]


results_dcc_d <- array(NA, dim = c(length(sel),20))
resultsCI_dcc_d <- array(NA, dim = c(length(sel),40))

## loop through sel participants
total <- length(sel)
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

for(i in 1:length(sel)){
  
#for(i in redo){
    partner <- unique(ddip.lag[ddip.lag$ID_DYAD==sel[i], "ID_INDIV"])
    ## Extract Time series per invididul in dyad
    r1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1],'pa_rel_std']
    r2 <- ddip.lag[ddip.lag$ID_INDIV==partner[2],'pa_rel_std']
    r <- scale( cbind(r1, r2) )
    fit1 = bmgarch(data = r, parameterization = "DCC", iterations = 1500, standardize_data = FALSE)
    
    pst <- rstan::summary(fit1$model_fit, c('c_h', 'a_h', 'b_h', 'a_q', 'b_q', 'S[1,2]', 'phi0', 'phi', 'theta', 'lp__'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')]
    
    ## Overwrite parameter if Rhat < 1.1
    results_dcc_d[i, ] <- pst[,1]*ifelse(pst[,'Rhat']>3, NA, 1)
    colnames(results_dcc_d) <- rownames(pst)
    ## Record lower and upper CI
    resultsCI_dcc_d[i,] <-  c(t(pst[,2:3]))
    low <- paste0(rownames(pst), '2.5%')
    hi <- paste0(rownames(pst), '97.5%')
    colnames(resultsCI_dcc_d)  <-   c(t(cbind(low, hi)))    
    ## Progress Bar
    setTxtProgressBar(pb, i)
}
    close(pb)
    
    ## focus on b2 for now?
    
  results_dcc_d
  resultsCI_dcc_d
  
  redo <- which(is.na(rowSums(results_dcc_d)))
  length(redo)
  redo
  total
  

save(results_dcc_d, file = 'results_dyad_dcc.RDat', compress = 'xz')
save(resultsCI_dcc_d, file = 'resultsCI_dyad_dcc.RDat', compress = 'xz')





## Predict relsatisfact
    sel
    names(ddip.lag)

    relsat <- NA
    reldur <- NA
    bup = NA
    for( i in 1:length(sel) ){
        relsat[i] <- mean(ddip.lag[ddip.lag$ID_DYAD==sel[i], 'relsatisfactc'], na.rm = T)
        reldur[i] <- mean(ddip.lag[ddip.lag$ID_DYAD==sel[i], 'totyrinrel'], na.rm = T)
        bup[i] = mean(ddip.lag[ddip.lag$ID_DYAD==sel[i], 'BREAKUP'], na.rm = T)
        
    }
    relsat

    bup[is.na(bup)] = 1
    bup
    
    dim(ddip.lag)
    dim(results)

    results_dcc_d

    summary( lm(relsat ~ results_bekk[,1:13]) )
    summary( lm(reldur ~ results_bekk[,1:13]) ) 
    summary( lm(relsat ~ results_dcc_d[,1:9]) )
    summary( lm(reldur ~ results_dcc_d[,1:9]) )

summary( glm(bup ~ results_bekk[,1:12], family = binomial) )
summary( glm(bup ~ results_dcc_d[,1:9], family = binomial) )

    
round(cor( cbind(relsat, results_bekk[,1:14]), use = "pairwise.complete.obs")[,1], 2)
round( cor( cbind(reldur, results_bekk), use = "pairwise.complete.obs")[,1] , 2 )
round( cor( cbind(relsat, results_dcc_d), use = "pairwise.complete.obs")[,1] , 2 ) 



    
#############################
## DCC
results_dcc <- array(NA, dim = c(length(sel),16))
resultsCI_dcc <- array(NA, dim = c(length(sel),32))

## loop through sel participants
total <- length(sel)
# create progress bar
pb <- tkProgressBar(title = "progress bar", min = 0, max = total, width = 300)
for(i in 1:length(sel)){
    i = 34
    partner <- unique(ddip.lag[ddip.lag$ID_DYAD==sel[i], "ID_INDIV"])
    ## Extract Time series per invididul in dyad
    r1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1],'pa_rel']
    r2 <- ddip.lag[ddip.lag$ID_INDIV==partner[2],'pa_rel']
    r <- cbind(r1, r2)
    ## Provide init sigma
    sigma1 <- var(r)   
    standat <- list(T = nrow(r), rts = r, sigma1 = sigma1, nt = ncol(r), ahead = ahead)
    dcc_fit <- sampling(dcc_mod, data = standat, verbose = TRUE, iter = 3000, control = list(adapt_delta = 0.95), init_r = 1, chains = 4)
    pst <- summary(dcc_fit, c('c_h', 'a_h', 'b_h', 'a_q', 'b_q', 'S[1,2]', 'b0', 'b1', 'b2', 'lp__'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')]
    ## Overwrite parameter if Rhat < 1.1
    results_dcc[i, ] <- pst[,1]*ifelse(pst[,'Rhat']>1.1, NA, 1)
    colnames(results_dcc) <- rownames(pst)
    ## Record lower and upper CI
    resultsCI_dcc[i,] <-  c(t(pst[,2:3]))
    low <- paste0(rownames(pst), '2.5%')
    hi <- paste0(rownames(pst), '97.5%')
    colnames(resultsCI_dcc)  <-   c(t(cbind(low, hi)))    
    ## Progress Bar
    setTkProgressBar(pb, i, label=paste( round(i/total*100, 0),  "% done"))
}

close(pb)

## focus on b2 for now?

results_dcc
resultsCI_dcc

save(results_dcc, file = 'results_dyad_dcc.RDat', compress = 'xz')
save(resultsCI_dcc, file = 'resultsCI_dyad_dcc.RDat', compress = 'xz')

#load(file = 'results_dyad.RDat')

## Create filter for results with posterior mass outside of zero
nonzero_dcc <- ifelse(sign(resultsCI_dcc[, seq(1, 32, by = 2)]) + sign(resultsCI_dcc[, seq(2, 32, by = 2)]) == 0, NA, 1)

## Return nonzero results:
round(results_dcc * nonzero_dcc, 2)

## END DCC
############################

#old.res <- results

names(ddip.lag)
head(ddip.lag)

df <- data.frame(results)
df


partner <- unique(ddip.lag[ddip.lag$ID_DYAD==sel[2], "ID_INDIV"])
    r1 <- ddip.lag[ddip.lag$ID_INDIV==partner[1],'pa_rel']

library(Hmisc)
r1l <- Hmisc::Lag(r1)
lm(r1 ~ r1l)

cor(df[, c('A.2.2.', 'B.2.2.', 'b1.2.')], use = 'pair')

mean(df$b1.2., na.rm = T)




## Per Dyad
a <- aggregate(ddip.lag$totyrinrel, by = list(ddip.lag$ID_DYAD), FUN = mean)
a
df$yearinrel <- a$x

## Per individual
## V2 Relsat
b <- aggregate(ddip.lag$V2RelSat, by = list(ddip.lag$ID_INDIV, ddip.lag$ID_DYAD), FUN = mean)
b
partner.1 <- array(NA, dim = c(length(sel), 1))
partner.2 <- partner.1
for(i in 1:length(sel)){
    partner.1[i] <- b[b$Group.2==sel[i],'x'][1]
    partner.2[i] <- b[b$Group.2==sel[i],'x'][2]
}
df$relsat.1 <- partner.1
df$relsat.2 <- partner.2

## V2 Anxiety (all NA)
##
names(ddip.lag)
b <- aggregate(ddip.lag$V1Anxiety, by = list(ddip.lag$ID_INDIV, ddip.lag$ID_DYAD), FUN = mean, na.action = na.omit)
b
partner.1 <- array(NA, dim = c(length(sel), 1))
partner.2 <- partner.1
for(i in 1:length(sel)){
    partner.1[i] <- b[b$Group.2==sel[i],'x'][1]
    partner.2[i] <- b[b$Group.2==sel[i],'x'][2]
}
df$anxiety.1 <- partner.1
df$anxiety.2 <- partner.2
df$anxiety.1

##
df[c(29,33),]

cor(na.omit(df[,c('Const.1.1.', 'Const.2.2.')]))

cor(na.omit(df[,c('corC.1.2.', 'A.1.2.', 'B.1.2.', 'b0.1.', 'b1.2.')]))

cor(na.omit(df[,c('A.1.2.', 'A.2.1.', 'B.1.2.', 'B.2.1.', 'b0.1.', 'b1.2.', 'relsat.1', 'relsat.2')]))

summary(lm(relsat.2 ~ A.1.2., data = df))

plot(df$b2.1.)

## Variance for women
mean(df$Const.1.1., na.rm = T)
## Variance for men
mean(df$Const.2.2., na.rm = T)
 cor(df)


cor(na.omit(df[,c('Const.1.1.', 'Const.2.2.','Const.1.2.', 'corC.1.2.',
          'A.1.1.', 'A.1.2.', 'A.2.1.', 'relsat.1', 'relsat.2')]))


df$anxiety.1

round(colMeans(df, na.rm = T), 2)

cor(na.omit(df[,c('Const.1.1.', 'Const.2.2.', 'A.1.1.', 'A.2.2.', 'A.1.2.', 'A.2.1.', 'B.1.1.', 'B.2.2.','B.1.2.', 'B.2.1.', 'b0.1.', 'b1.2.', 'yearinrel', 'anxiety.1', 'anxiety.2')]))[, 14:15]

df

round(cov2cor(var(df, na.rm = TRUE )), 2)

summary(lm(A.2.2. ~ relsat.1 + relsat.2, data = df))


g <- array(c(1,5), dim = c(2,1))

cor(na.omit(cbind(abs(df$A.1.2. - df$A.2.1.),    df[c('relsat.1', 'relsat.2')])))
cor(na.omit(cbind(abs(df$B.1.2. - df$B.2.1.),    df[c('relsat.1', 'relsat.2')])))
