source('settings.r')

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
approximations = list(mra = mra, low.rank = low.rank)#, exact = exact)


XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)


##### scalable FFBS
# run smoothing for each sample; maybe this can be optimized
#samples = foreach( sample.no = 1:Nsamples ) %dopar% {

smoothingResults = FFBS('mra', XY$y)
samples = smoothingResults$samples
results = smoothingResults$filteringResults

meansPS = getMeansFromSamples(samples)
sdsPS   = getSDFromSamples(samples)
meansF  = getMeansFromFilter(results$preds)
sdsF    = getSDFromFilter(results$preds)
plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep = ""))


###### compare with naive Kalman Smoother
#KSsmooth = KalmanSmoother(XY$y)
#meansKS = getMeansFromSmoother( KSsmooth$smoothed )
#sdsKS   = getSDFromSmoother( KSsmooth$smoothed )
#plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, '-ps-ks', sep = ""))
