source("../advection.r")
source('../simulate-data.r')
source("../smoother.r")
source("../FFBS.r")
source("../filter.r")
source("sig2-sample.r")
library(RandomFields)
library(invgamma)
library(Matrix)

######### set parameters #########
#set.seed(1988)
N = 34**2
M = 50
NSAMPLES = 100
FRAC_OBS = 0.3
TMAX = 20
# 34 x 34 parameters
evolFun = function(X) evolDiff(X, adv = 0.01, diff = 0.00004)
#evolFun = function(X) evolDiff(X, adv = 0.1, diff = 0.0004)


## covariance parameters
SIG_02 = 1
SIG2 = 0.1; RANGE = .15; SMOOTH = 0.5; 
COVPARMS = c(SIG2, RANGE, SMOOTH)
PRIOR_COVPARMS = c(SIG_02, RANGE, SMOOTH)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs), COVPARMS)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, PRIOR_COVPARMS)

## likelihood settings
ME_VAR = 0.05;
DATA_MODEL = "gauss"  
LIK_PARMS = list(data.model = DATA_MODEL, sigma = sqrt(ME_VAR))
ALPHA = 0.001
BETA = 0.001
