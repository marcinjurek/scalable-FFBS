source('/home/marcin/MRS/gibbs-sampler/settings.r')

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
Sig0Model = RandomFields::RMmatern(nu = smooth, scale = range)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
Sig0Mat = RandomFields::RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)


RandomFields::RFoptions(storing = FALSE)
x0 = matrix( x0, ncol = 1)
zeros = matrix(rep(0, n), ncol = 1)


## define Vecchia approximation
#exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
#low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra)#, low.rank = low.rank)#, exact = exact)

Qinv = solve(Sig0Mat)
XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

###### naive Kalman Smoother
#KSsmooth = KalmanSmoother(XY$y)

##### scalable FFBS
# run smoothing for each sample; maybe this can be optimized
#samples = foreach( sample.no = 1:Nsamples ) %dopar% {
#scores = foreach(iter = 1:max.iter) %dopar% {
#scores = list()

sig2draw = 0.5*sig2#rinvgamma(1, alpha, beta)
cat(paste("Sig2 =", sig2draw, '\n'))

for (iter in 1:max.iter) {
    cat(paste("Starting iteration ", iter, "\n", sep = ""))
    local.scores = list()    
    for (name in names(approximations)) {
        cat(paste("Using approximation ", name, "\n", sep = ""))
        
        draws = list(state = list(), sig2 = rep(0, Nsamples))
        
        for (sample.no in 1:Nsamples) {
            cat(paste("Sample no", sample.no, "\n"))
            smoothingResults = FFBS(name, XY$y, Nsamples = 1, sig2draw)
            sampledStates = smoothingResults$samples[[1]]
            
            sig2draw = sig2.sample(sampledStates, evolFun, alpha, beta, Qinv)
            cat(paste("Sig2 =", sig2draw, '\n'))
            
            draws$state[[ sample.no ]] = sampledStates
            draws$sig2[ sample.no ] = sig2draw
        }

    }
   # scores[[iter]] = local.scores
}

pdf("trajectory-sig2.pdf")
plot(draws$sig2)
dev.off()


#scores.KS = getScores(XY$x, smoothed = KSsmooth$smoothed)
# scores.per.method = list()
# for (name in names(approximations)) {
#     scores.per.method[[ name ]] = rep(0, Tmax)
#     for (iter in 1:max.iter) {
#         scores.per.method[[ name ]] = scores.per.method[[ name ]] + scores[[ iter ]][[ name ]]/max.iter
#     }
# }
# 
# pdf("CRPS.pdf")
# line.colours = list(mra = "red", exact = "black", low.rank = "blue")
# plot(NULL, xlim = c(1, Tmax), ylim = range(sapply(scores, range)), ylab = "CRPS", xlab = "time")
# for (name in names(approximations)) {
#     lines(scores.per.method[[name]], col = line.colours[[name]])
# }
# legend("topright", names(approximations), col = c("red", "black", "blue"), lwd = c(1, 1, 1))
# dev.off()

#meansPS = getMeansFromSamples(samples)
#sdsPS   = getSDFromSamples(samples)
#meansF  = getMeansFromFilter(results$preds)
#sdsF    = getSDFromFilter(results$preds)
#plot2Results(XY, meansPS, sdsPS, meansF, sdsF, paste(data.model, '-ps-f', sep = ""))



#meansKS = getMeansFromSmoother( KSsmooth$smoothed )
#sdsKS   = getSDFromSmoother( KSsmooth$smoothed )
#plot2Results(XY, meansPS, sdsPS, meansKS, sdsKS, paste(data.model, '-ps-ks', sep = ""))
