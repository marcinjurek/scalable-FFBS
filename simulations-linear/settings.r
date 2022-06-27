source("../aux-functions.r")
source("../advection.r")
source('../simulate-data.r')
source("../smoother.r")
source("../FFBS.r")
source("../filter.r")
library(doParallel)
library(invgamma)
library(RandomFields)
library(Matrix)
resultsDir = "simulations-linear"
NCORES = 5
cl = makeCluster(NCORES, type="FORK")  
registerDoParallel(cl) 

######### set parameters #########
set.seed(1988)
spatial.dim = 2
#n = 20 ** 2
#m = 30
n = 34**2
m = 50
Nsamples = 50
frac.obs = 0.3
Tmax = 20
max.iter = 50
diffusion = 0.00004
advection = 0.01
#evolFun = function(X) X
evolFun = function(X) evolDiff(X, adv = advection, diff = diffusion)
options(digits=4)

## covariance parameters
sig_02 = 1
sig2 = 0.1; range = .15; smooth = 0.5; 
covparms = c(sig2, range, smooth)
prior_covparams = c(sig_02, range, smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparams)

me.var = 0.05
data.model = "gauss"
lik_params = list(data.model = data.model, sigma = sqrt(me.var), alpha = 2)
