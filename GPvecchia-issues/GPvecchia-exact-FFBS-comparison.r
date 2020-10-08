### This script shows that using the FFBS approach and GPvecchia
### does not produce correct estimates of the variance
###
### Below, within the for loop we simulate the field of interest, x,
### and the observation vector y = Hx + v, where v ~ N(0, me.var).
### We set some very small measurement error so that our data
### is very close to the true field. The H matrix simply selects certain
### elements from x (specifically, it selects frac.obs=50% of them)
###
### Next we generate a synthetic field x.hat and synthetic data y.hat
### in the same way. x.hat and x have the same distirbution as do
### y.hat and y. Then, after specifying the vecchia approximation,
### we use y-y.hat as observations and use the GPvecchia function to
### calculate the posterior distribution.
###
### As the last step we add the synthetic state vector x.hat back in.
### This entire procedure corresponds to the FFBS algorithm for a model
### with Tmax=1.
###
### At the end of each iteration we calculate the estimated variance
### parameter. We expect the mean of these estimates to be similar to the
### true value sig2=0.5. This is the case when we use the exact posterior
### distirbution (set approx=exact), but it is not the case when approx=mra.
###
### For Niter = 200, mean(est.var) = 0.5510268 and sd(est.var) = 0.0244143
### when approx = mra. On the contrary, when approx = exact with just Niter=100
### we have mean(est.var) = 0.04994163 and sd(est.var) = 0.02050539.
### A simple t-test shows that the sampling error is unlikely to be behind
### the elevated value of the MRA-based estimate.

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

Niter = 200 # number of repeti

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

    if(iter %% 10==0){
        cat(paste("iter", iter, "\n"))
    }

    ### truth
    x = Matrix( RFsimulate(model = Sig0, x = locs[,1], y = locs[,2], spConform = FALSE) )
    v = Matrix( me.var*rnorm(n) )
    
    y = H %*% (x + v)

    ### synthetic
    x.hat = Matrix( RFsimulate(model = Sig0, x = locs[,1], y = locs[,2], spConform = FALSE) )
    v.hat = Matrix( me.var*rnorm(n) )
    y.hat = H %*% (x.hat + v.hat)
    

    ### calculate stats

    # the simple way
    #mu.post.simple = CovXY %*% PrecY %*% (y - y.hat)

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
    obs.aux[obs.inds] = y - y.hat

    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = rep(0, n),
                                       likelihood_model = data.model,
                                       covmodel = covmodel, covparms = c(sig2, range, smooth),
                                       likparms = lik.params, return_all = TRUE)
    mu.post = matrix(preds.aux$mean, ncol = 1)
    # ----------------
    
    # generate sample
    x.sample = mu.post + x.hat

    est.mu[iter] = as.numeric(mean(x.sample))
    est.var[iter] = as.numeric(t(x.sample) %*% (PrecX*sig2) %*% x.sample / (n-1))
}

#hist(est.mu)
hist(est.var)  
