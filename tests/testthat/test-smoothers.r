#Source('../aux-functions.r')
source('../../simulate-data.r')
source("../../smoother.r")
source("../../FFBS.r")
source("../../filter.r")
library(RandomFields)
library(Matrix)

######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 20**2
m = 30
maxiter = 1
frac.obs = 0.3
Tmax = 2
evolFun = function(X) X

## covariance parameters
sig_02 = 1
sig2 = 0.25; range = .25; smooth = 0.5; 
covparms = c(sig2,range,smooth)
prior_covparms = c(sig_02, range, smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparms)

## likelihood settings
me.var = 0.1;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var))


######### initialize #########
## generate grid of pred.locs
grid.oneside = seq(0, 1, length = round(sqrt(n)))
grid = expand.grid(grid.oneside, grid.oneside)
locs = as.matrix(grid) 

## set initial state
Sig0Model = RMwhittle(nu = smooth, scale = range, var = 1)
x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))
RFoptions(storing = FALSE)
Sig0Mat = RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)

## define Vecchia aproximation
#mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = "firstm")
approximations = list(exact = exact)#, mra = mra)



XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                 sig2 = sig2, smooth = smooth, range = range, locs = locs)

app = approximations[["exact"]]
smoothed = smoothedMeans(app, XY$y, lik.params, prior_covparams, covparms, evolFun)
Ksmoothed = KalmanSmoother(XY$y, lik.params, Sig0Mat, sig2*Sig0Mat, evolFun)
    
diffs = list()
diffsf = list()
for (t in 1:Tmax) {
    diffs[[t]] = sum(abs(smoothed$means[[t]] - Ksmoothed$smoothed[[t]]))
    print(diffs[[t]])
}
maxdiff = max(unlist(diffs))
cat(sprintf("Maximum difference between the exact vecchia filter and the Kalman smoother is %f\n", maxdiff))
