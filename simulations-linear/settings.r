setwd("~/MRS/simulations-linear")
source('aux-functions.r')


######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 34**2
m = 50
frac.obs = 0.1
Tmax = 5
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Nsamples = 10


## covariance parameters
sig2 = 0.5; range = .1; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.05;
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
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)
