source("../../advection.r")
source('../../simulate-data.r')
source("../../smoother.r")
source("../../FFBS.r")
source("../../filter.r")
library(RandomFields)
library(Matrix)

######### set parameters #########
#set.seed(1988)
spatial.dim = 2
N = 34**2
m = 50
maxiter = 5
frac.obs = 0.3
Tmax = 20
evolFun = function(X) evolDiff(X, adv = 0.01, diff = 0.00004)
#evolFun = function(X) evolDiff(X, adv = 0.1, diff = 0.0004)


## covariance parameters
sig_02 = 1
sig2 = 0.1; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
prior_covparms = c(sig_02, range, smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparms)

## likelihood settings
me.var = 0.05;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var))
