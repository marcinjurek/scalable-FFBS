setwd("~/FFBS/gibbs-sampler")
source('../aux-functions.r')
source('../smoother.r')
source('../simulate-data.r')
source('../advection.r')
source('../filter.r')
source('../FFBS.r')
source('sig2-sample.r')
resultsDir = "gibbs-sampler"
library(RandomFields)
library(Matrix)
library(invgamma)

options(digits=4)

######### set parameters #########
#set.seed(1988)
spatial.dim = 2
n = 20**2
m = 30
frac.obs = 1.0
Tmax = 2
advection = 0
diffusion = 0


evolFun = function(X) X#function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Nsamples = 100
max.iter = 1


## covariance parameters
sig_02 = 1
sig2 = 0.25; range = .25; smooth = 0.5; 
covparms = c(sig2, range, smooth)
prior_covparams = c(sig_02, range, smooth)

covfun.base = function(locs, covparms) GPvecchia::MaternFun(fields::rdist(locs), covparms)
covfun.base.d = function(D, covparms) GPvecchia::MaternFun(D, covparms)

covfun = function(locs) covfun.base(locs, covparms)
covfun.d = function(D) covfun.base.d(D, covparms)


## likelihood settings
me.var = 0.1;
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
alpha = 5
beta = 1

mu.a = advection
sig2.a = mu.a/6
