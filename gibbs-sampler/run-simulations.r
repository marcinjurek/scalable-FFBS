source("settings.r")

######### initialize #########
## generate grid of pred.locs
grid.oneside = seq(0, 1, length = round(sqrt(N)))
grid = expand.grid(grid.oneside, grid.oneside)
locs = as.matrix(grid) 

## set initial state
Sig0Model = RMwhittle(nu = SMOOTH, scale = RANGE, var = SIG_02)
x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))
RFoptions(storing = FALSE)
Sig0Mat = RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)

## define Vecchia aproximation
mra = GPvecchia::vecchia_specify(locs, M, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, N - 1, conditioning = "firstm")
lowrank = GPvecchia::vecchia_specify(locs, M, conditioning = 'firstm')
#approximations = list(lowrank = lowrank, exact = exact, mra = mra)
approximations = list(exact = exact)


Linvs = list(truth = solve(t(chol(Sig0Mat))))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}


XY = simulate.xy(x0, evolFun, NULL, FRAC_OBS, LIK_PARMS, TMAX,
                 sig2 = SIG2, smooth = SMOOTH, range = RANGE, locs = locs)

L = Linvs[["truth"]]
sse = 0
for (t in 2:TMAX) {
    xs1 = XY$x[[t - 1]]
    xs2 = XY$x[[t]]
    w = xs2 - evolFun(xs1)
    sse = sse + t(w) %*% t(L) %*% L %*% w
}
est = as.numeric(sse / ((N - 1) * (TMAX - 1)))
cat(sprintf("%s Estimate from true field %f\n", Sys.time(), est))


sigmasL = list()
for (name in names(approximations)) {

    cat(sprintf("%s Sampling using %s\n", Sys.time(), name))
    app = approximations[[name]]
    sig2_1 = runif(1)
    sigmas = rep(sig2_1, NSAMPLES)
    covparams = COVPARMS

    for (sample_no in 2:NSAMPLES) {

        if (sample_no %% 10 == 0) {
            cat(sprintf("%s\t working on sample %d\n", Sys.time(), sample_no))
        }
        
        smoothed = FFBS(app, XY$y, LIK_PARMS, PRIOR_COVPARMS, covparams, evolFun)
        states = smoothed$samples[[1]]
        L = Linvs[[name]]
        sigmas[sample_no] = sig2_sample(states, evolFun, ALPHA, BETA, Linvs[[name]])
        covparams[1] = sigmas[sample_no]
    }

    sigmasL[[name]] = sigmas
}

sigmasL = as.data.frame(sigmasL)
write.csv(sigmasL, file = "gibbs-samples.csv")
plot(sigmasL$exact, ylim = range(sigmasL), type = "l");
#lines(sigmasL$mra, col = "red");
#lines(sigmasL$lowrank, col = "blue")
