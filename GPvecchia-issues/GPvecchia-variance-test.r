rm(list=ls())

source('simulations-linear/getMatCov.r')
Rcpp::sourceCpp('simulations-linear/getMatCovFromFactor.cpp')
suppressPackageStartupMessages(library(fields))
suppressPackageStartupMessages(library(RandomFields, quietly = T))
suppressPackageStartupMessages(library(Matrix))


### settings ###
set.seed(4796)
sqrt.n = 34
n = sqrt.n**2

m = 50
sig2 = 0.5
smooth = 0.5
range = 0.2
me.var = 1e-8
frac.obs = 0.5

Niter = 100 # number of repeti

### initialize fields ###
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 

Sig0 = RMmatern(nu = smooth, scale = range, var = sig2)

obs.inds = sample(1:n, round(frac.obs*n))
H = Matrix::Diagonal(n)[obs.inds,]

### plot fields ###
#quilt.plot(locs, as.numeric(y), nx=sqrt.n, ny=sqrt.n, main='data')
#quilt.plot(locs, as.numeric(x), nx=sqrt.n, ny=sqrt.n, main='truth')

### precompute matrices for FFBS ###
CovX = Matrix( RFcovmatrix(Sig0, locs) )
#CovY = H %*% (CovX + me.var*Matrix::Diagonal(n)) %*% t(H)
#CovXY = CovX %*% t(H)
PrecX = solve(CovX)
#PrecY = solve(CovY)

est.var = est.mu = rep(NA, Niter)

for(iter in 1:Niter){


    ### truth
    x = Matrix( RFsimulate(model = Sig0, x = locs[,1], y = locs[,2], spConform = FALSE) )
    v = Matrix( me.var*rnorm(n) )
    
    y = H %*% (x + v)

    
    # using GPvecchia
    # ----------------
    data.model = 'gauss'   
    lik.params = list(data.model = data.model, sigma = sqrt(me.var))

    #exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose=FALSE)
    mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose=FALSE)
    approx = mra
    
    covfun.d = function(D) GPvecchia::MaternFun(D, c(sig2, range, smooth))
    #covfun = function(locs1, locs2) GPvecchia::MaternFun(rdist(locs1, locs2), c(sig2, range, smooth))
    covmodel = getMatCov(approx, covfun.d)

    obs.aux = rep(NA, n)
    obs.aux[obs.inds] = y

    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = rep(0, n),
                                       likelihood_model = data.model,
                                       covmodel = covmodel, covparms = c(sig2, range, smooth),
                                       likparms = lik.params, return_all = TRUE)
    mu.post = matrix(preds.aux$mean, ncol = 1)
    # ----------------
    
    # generate sample
    est.mu[iter] = as.numeric(mean(mu.post))
    est.var[iter] = as.numeric(t(mu.post) %*% (PrecX*sig2) %*% mu.post / (n-1))
}

print(summary(est.mu))
print(summary(est.var))
