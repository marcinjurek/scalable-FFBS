setwd("~/FFBS/simulations-linear")
source('aux-functions.r')
source('simulate-data.r')
source('smoother.r')
source('filter.r')
source('advection.r')
Rcpp::sourceCpp('getMatCovFromFactor.cpp')
source('getMatCov.r')
source('FFBS.r')
library(Matrix)


######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 34**2
m = 30
frac.obs = 0.3
Tmax = 2
#diffusion = 0.00004
#advection = 0.1
evolFun = function(X) X
#evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Nsamples = 20

## covariance parameters
sig_02 = 1
sig2 = 0.25; range = .1; smooth = 0.5; 
covparms = c(sig2,range,smooth)
prior_covparams = c(sig_02, range, smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 1e-16;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var))
