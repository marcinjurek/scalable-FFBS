setwd("~/MRS/gibbs-sampler")
source("results-display.r")
source('aux-functions.r')
source('smoother.r')
source('getMatCov.r')
source('simulate-data.r')
source('advection.r')
source('FFBS.r')
source('sig2-sample.r')
source('a-sample.r')
Rcpp::sourceCpp('getMatCovFromFactor.cpp')
resultsDir = "gibbs-sampler"
library(Matrix)
library(doParallel)
library(invgamma)
registerDoParallel(cores = 6)


options(digits=4)


######### set parameters #########
#set.seed(1988)
spatial.dim = 2
n = 34**2
m = 34
frac.obs = 0.5
Tmax = 1
diffusion = 0.0#0004
advection = 0.0#1


evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Nsamples = 1
max.iter = 10


## covariance parameters
sig2 = 0.5; range = 0.5; smooth = 0.5; 
covparms = c(sig2,range,smooth)

covfun.base = function(locs, covparms) GPvecchia::MaternFun(fields::rdist(locs), covparms)
covfun.base.d = function(D, covparms) GPvecchia::MaternFun(D, covparms)

covfun = function(locs) covfun.base(locs, covparms)
covfun.d = function(D) covfun.base.d(D, covparms)


## likelihood settings
me.var = 1e-8;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"  
}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha = 2)

## parameter prior settings
alpha = 0.1
beta = 0.1
