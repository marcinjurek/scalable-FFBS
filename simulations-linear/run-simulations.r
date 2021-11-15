source('settings.r')

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
Sig0Model = RandomFields::RMwhittle(nu = smooth, scale = range, var = sig_02)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
RandomFields::RFoptions(storing = FALSE)
x0 = matrix(x0, ncol = 1)

## define Vecchia approximation
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra, low.rank = low.rank, exact = exact)

PrecMat = list()
for (name in names(approximations)) {
    appr = approximations[[name]]
    L =  getL00(appr, covfun.d, locs) / sqrt(sig2)
    PrecMat[[name]] = solve(L %*% t(L))
}

XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

##### scalable FFBS
# run smoothing for each sample; maybe this can be optimized
#amples = foreach( sample.no = 1:Nsamples ) %dopar% {

for (name in names(approximations)) {
    print(name)
    smoothingResults = FFBS(appr, XY$y, lik.params, prior_covparams)
    samples = smoothingResults$samples

    vars = rep(0, Nsamples)
    for (sampleNo in 1:Nsamples) {
        xf1 = samples[[sampleNo]][[1]]
        xf2 = samples[[sampleNo]][[2]]
        ef = xf2 - xf1
        est_var = t(ef) %*% PrecMat[["mra"]] %*% ef / (n - 1)
        vars[sampleNo] = est_var
    }
    print(vars)
}
        #results = smoothingResults$filteringResults

        #meansPS = getMeansFromSamples(samples)
        #sdsPS   = getSDFromSamples(samples)
        #meansF  = getMeansFromFilter(results$preds)
        #sdsF    = getSDFromFilter(results$preds)
        #plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep = ""))


###### compare with naive Kalman Smoother
#KSsmooth = KalmanSmoother(XY$y)
#meansKS = getMeansFromSmoother( KSsmooth$smoothed )
#sdsKS   = getSDFromSmoother( KSsmooth$smoothed )
#plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, '-ps-ks', sep = ""))
