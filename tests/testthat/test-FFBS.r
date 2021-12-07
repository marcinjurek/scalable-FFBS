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
maxiter = 100
frac.obs = 1.0
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
me.var = 1e-10;
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
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = "firstm")
approximations = list(exact = exact, mra = mra)


Linvs = list(truth = solve(chol(Sig0Mat)))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}



d = data.frame(truth = c(), exact = c(), mra = c())
for (iter in 1:maxiter) {
    cat(sprintf("======== iter %d ========\n", iter))
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                     sig2 = sig2, smooth = smooth, range = range, locs = locs)
    es = list()
    es[["truth"]] = XY$x[[2]] - XY$x[[1]]

    for (name in names(approximations)) {
        app = approximations[[name]]
        filtered = VecchiaFilter(app, XY$y, lik.params, prior_covparms, covparms, evolFun)
        smoothed = FFBS(app, XY$y, lik.params, prior_covparams, covparms, evolFun, Num_samples = 1, covparams = NULL)
        xs1 = smoothed$samples[[1]][[1]]
        xs2 = smoothed$samples[[1]][[2]]
        #xs1 = filtered$preds[[1]]
        #xs2 = filtered$preds[[2]]
        es[[name]] = xs2 - xs1
    }

    est = list()
    for (name in names(es)) {
        e = es[[name]]
        L = Linvs[[name]]
        est[[name]] = as.numeric(t(e) %*% t(L) %*% L %*% e / (n - 1))
    }
    dloc = data.frame(truth = est[["truth"]], exact = est[["exact"]], mra = est[["mra"]])
    d = rbind(d, dloc)
}

oldpar = par(mfrow = c(3, 1))
hist(d$truth, xlim = range(d), breaks = 10)
hist(d$exact, xlim = range(d), breaks = 10)
hist(d$mra, xlim = range(d), breaks = 10)
par(oldpar)
