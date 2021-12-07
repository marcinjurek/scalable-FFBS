#source('../aux-functions.r')
source('~/FFBS/simulate-data.r')
source("~/FFBS/filter.r")
library(RandomFields)
library(Matrix)

######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 3**2
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

## define Vecchia approximation
vecchia = GPvecchia::vecchia_specify(locs, n - 1, conditioning = "firstm")



XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                 sig2 = sig2, smooth = smooth, range = range, locs = locs)


Q_cov = sig2 * Sig0Mat
Vfiltered = VecchiaFilter(vecchia, XY$y, lik.params, prior_covparms, covparms, evolFun)
Kfiltered = KalmanFilter(XY$y, lik.params, Sig0Mat, Q_cov, evolFun)

Kp = Kfiltered$preds
Vp = Vfiltered$preds
VLs = Vfiltered$Ls
VVs = Vfiltered$Vs

#VL = Vfiltered$Ls[[1]]
#KL = Kfiltered$Ls[[1]]




#oldpar = par(mfrow = c(1, 3))
#fields::quilt.plot(locs[, 1], locs[, 2], XY$x[[1]], nx = 20, ny = 20, main = "truth")
#fields::quilt.plot(locs[, 1], locs[, 2], Kfiltered$preds[[1]], nx = 20, ny = 20, main = "Kalman")
#fields::quilt.plot(locs[, 1], locs[, 2], Vfiltered$preds[[1]], nx = 20, ny = 20, main = "Vecchia exact")
#par(oldpar)
