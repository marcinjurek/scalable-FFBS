setwd("~/FFBS/simulations-linear")
source('../aux-functions.r')
source('../simulate-data.r')
#Rcpp::sourceCpp('getMatCovFromFactor.cpp')
#source('getMatCov.r')
source("../filter.r")
library(Matrix)
library(RandomFields)

######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 20**2
m = 30
frac.obs = 1.0
Tmax = 2
evolFun = function(X) X

## covariance parameters
sig2 = 0.25; range = .25; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.1;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var))
