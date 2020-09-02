setwd("~/MRS/simulations-linear")
source("results-display.r")
source('aux-functions.r')
source('smoother.r')
source('getMatCov.r')
source('settings.r')
source('simulate-data.r')
source('advection.r')
Rcpp::sourceCpp('getMatCovFromFactor.cpp')
resultsDir = "simulations-linear"
library(Matrix)
library(doParallel)
registerDoParallel(cores = 6)



## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
Sig0Model = RandomFields::RMmatern(nu = smooth, scale = range)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
RandomFields::RFoptions(storing = FALSE)
x0 = matrix( sig2 * x0, ncol = 1)
zeros = matrix(rep(0, n), ncol = 1)


## define Vecchia approximation
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)



# run smoothing for each sample; maybe this can be optimized

#scores = foreach( iter=1:max.iter) %dopar% {
#samples = foreach( sample.no = 1:Nsamples ) %dopar% {
samples = list()
for (sample.no in 1:Nsamples) {

    cat(paste("Working on sample no", sample.no, "\n"))
    
    XYplus = simulate.xy(zeros, evolFun, NULL, XY$y, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    Y = mapply('-', XY$y, XYplus$y, SIMPLIFY = FALSE)
    results = smoothedMeans('mra', Y)
    
    smeans = results[[ "means" ]]
    sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)
    cat("done working on sample\n")
    samples[[sample.no]] = sample
}



meansPS = getMeansFromSamples(samples)
sdsPS   = getSDFromSamples(samples)
meansF  = getMeansFromFilter(results$preds)
sdsF    = getSDFromFilter(results$preds)
plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep=""))


KFsmooth = KalmanSmoother(XY$y)
meansKS = getMeansFromSmoother( KFsmooth$smoothed )
sdsKS   = getSDFromSmoother( KFsmooth$smoothed )
plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, 'gauss-ps-ks', sep=""))
