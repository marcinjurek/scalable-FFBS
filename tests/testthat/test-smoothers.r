source("test-settings.r")

######### initialize #########
## generate grid of pred.locs
grid.oneside = seq(0, 1, length = round(sqrt(N)))
grid = expand.grid(grid.oneside, grid.oneside)
locs = as.matrix(grid) 

## set initial state
Sig0Model = RMwhittle(nu = smooth, scale = range, var = 1)
x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))
RFoptions(storing = FALSE)
Sig0Mat = RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)

## define Vecchia aproximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, N - 1, conditioning = "firstm")
approximations = list(exact = exact, mra = mra)


XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                 sig2 = sig2, smooth = smooth, range = range, locs = locs)


smoothed = list()
for (name in names(approximations)) {
    app = approximations[[name]]
    smoothed[[name]] = smoothedMeans(app, XY$y, lik.params, prior_covparms, covparms, evolFun)
}
Ksmoothed = KalmanSmoother(XY$y, lik.params, Sig0Mat, sig2*Sig0Mat, evolFun)


oldpar = par(mfrow = c(2, 2))
fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x[[Tmax]]), nx = 20, ny = 20, main = "truth")
fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(Ksmoothed$smoothed[[Tmax]]), nx = 20, ny = 20, main = "Kalman")
fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(smoothed[["exact"]][["means"]][[Tmax]]), nx = 20, ny = 20, main = "Vecchia exact")
fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(smoothed[["mra"]][["means"]][[Tmax]]), nx = 20, ny = 20, main = "MRA")
par(oldpar)



Linvs = list(truth = solve(t(chol(Sig0Mat))))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}


sse = list()
sse[["truth"]] = 0
sse[["kalman"]] = 0
Qinv = Linvs[["truth"]]
for (t in 2:Tmax) {
    w = XY$x[[t]] - evolFun(XY$x[[t - 1]])
    wK = Ksmoothed$smoothed[[t]] - evolFun(Ksmoothed$smoothed[[t - 1]])
    sse[["truth"]]  = sse[["truth"]]  + (t(w) %*% Matrix::t(Qinv) %*% Qinv %*% w)
    sse[["kalman"]] = sse[["kalman"]] + (t(wK) %*% Matrix::t(Qinv) %*% Qinv %*% wK)
}
sse[["truth"]] = as.numeric(sse[["truth"]] / ((N - 1)*(Tmax - 1)))
sse[["kalman"]] = as.numeric(sse[["kalman"]] / ((N - 1)*(Tmax - 1)))


for (name in names(approximations)) {
    sse[[name]] = 0
    Qinv = Linvs[[name]]
    states = smoothed[[name]][["means"]]
    for (t in 2:Tmax) {
        w = states[[t]] - evolFun(states[[t - 1]])
        sse[[name]] = sse[[name]] + (t(w) %*% Matrix::t(Qinv) %*% Qinv %*% w)
    }
    sse[[name]] = as.numeric(sse[[name]]/((N - 1)*(Tmax - 1)))
}
print(sse)


diffs = list()
diffsf = list()
for (t in 1:Tmax) {
    diffs[[t]] = sum(abs(smoothed[["exact"]][["means"]][[t]] - Ksmoothed$smoothed[[t]]))
}
maxdiff = max(unlist(diffs))
cat(sprintf("Maximum difference between the exact vecchia smoother and the Kalman smoother is %f\n", maxdiff))
