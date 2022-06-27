cat(paste(Sys.time(), " Calculating CRPS started\n", sep = ""))
source('settings.r')
library(RandomFields)
library(doParallel)
registerDoParallel(cores = NCORES)


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
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = "firstm")
low_rank = GPvecchia::vecchia_specify(locs, m, conditioning = 'firstm')


approximations = list(mra = mra, lr = low_rank, exact = exact)



###### naive Kalman Smoother
#KSsmooth = KalmanSmoother(XY$y)

##### scalable FFBS
# run smoothing for each sample; maybe this can be optimized
#samples = foreach( sample.no = 1:Nsamples ) %dopar% {
#scores = foreach(iter = 1:max.iter) %dopar% {
scores = list()
for (iter in 1:max.iter) {
    cat(paste(Sys.time(), " Starting iteration ", iter, "\n", sep = ""))
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik_params, Tmax,
                     sig2 = sig2, smooth = smooth, range = range, locs = locs)
    local_scores = list()    
    for (name in names(approximations)) {
        cat(paste(Sys.time(), " Using approximation ", name, "\n", sep = ""))
        smoothed = FFBS(approximations[[name]], XY$y, lik_params, prior_covparams,
                        covparms, evolFun, Num_samples = Nsamples)
        smoothedSamples = smoothed$samples
        local_scores[[name]] = getScores(XY$x, smoothedSamples)
    }
    scores[[iter]] = local_scores
    tosave = Reduce('+', lapply(scores, as.data.frame))
    write.csv(tosave, file = sprintf("sum-of-scores-after-%d-iterations.csv", iter))
}




#scores.KS = getScores(XY$x, smoothed = KSsmooth$smoothed)
#scores.per.method = list()
#scores.per.method[["exact"]] = rep(0, Tmax)
#for (iter in 1:max.iter) {
#    scores.per.method[["exact"]] = scores.per.method[["exact"]] + scores[[ iter ]][["exact"]]/max.iter
#}

#for (name in names(approximations)) {
#    if (name!="exact") {
#        scores.per.method[[ name ]] = rep(0, Tmax)
#        for (iter in 1:max.iter) {
#            scores.per.method[[ name ]] = scores.per.method[[ name ]] + scores[[ iter ]][[ name ]]/max.iter
#        }
#    }
#    scores.per.method[[name]] = scores.per.method[[name]] / scores.per.method[["exact"]]
#}

#pdf("CRPS.pdf")
#line.colours = list(exact = "black", mra = "red", lr = "blue")
#plot(NULL, xlim = c(1, Tmax), ylim = range(sapply(scores.per.method, range)), ylab = "CRPS", xlab = "time")
#for (name in names(approximations)) {
#    lines(scores.per.method[[name]], col = line.colours[[name]])
#}
#legend("topright", names(approximations), col = unlist(line.colours[names(approximations)]), lwd = c(1, 1, 1))
#dev.off()

#meansPS = getMeansFromSamples(samples)
#sdsPS   = getSDFromSamples(samples)
#meansF  = getMeansFromFilter(results$preds)
#sdsF    = getSDFromFilter(results$preds)
#plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep = ""))



#meansKS = getMeansFromSmoother( KSsmooth$smoothed )
#sdsKS   = getSDFromSmoother( KSsmooth$smoothed )
#plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, '-ps-ks', sep = ""))
