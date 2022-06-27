cat(sprintf("%s simulations started\n", Sys.time()))
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
#approximations = list(mra = mra)#, low.rank = low.rank, exact = exact)
approximations = list(mra = mra)#, low.rank = low.rank, exact = exact)



Linvs = list()#truth = solve(chol(Sig0Mat)))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}



XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik_params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)


for (t in 1:Tmax) {
    pdf(sprintf("simualted-field-%d.pdf", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x[[t]]), nx = sqrt(n), ny = sqrt(n))
    dev.off()
}


##### scalable FFBS
# run smoothing for each sample; maybe this can be optimized
#for( sample.no in 1:Nsamples ) {
#    #samples = foreach( sample.no = 1:Nsamples ) %dopar% {
#    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
#    smoothingResults = FFBS(mra, XY$y, lik.params, prior_covparams, covparms, evolFun, Num_samples = 1, covparams = NULL)
#    smoothingResults$samples
#
#    
#     #results = smoothingResults$filteringResults
#
#     #meansPS = getMeansFromSamples(samples)
#     #sdsPS   = getSDFromSamples(samples)
#     #meansF  = getMeansFromFilter(results$preds)
#     #sdsF    = getSDFromFilter(results$preds)
#     #plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep = ""))
# }

###### compare with naive Kalman Smoother
#KSsmooth = KalmanSmoother(XY$y)
#meansKS = getMeansFromSmoother( KSsmooth$smoothed )
#sdsKS   = getSDFromSmoother( KSsmooth$smoothed )
#plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, '-ps-ks', sep = ""))
