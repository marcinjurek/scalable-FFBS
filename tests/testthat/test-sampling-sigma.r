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
lowrank = GPvecchia::vecchia_specify(locs, m, conditioning = 'firstm')
approximations = list(exact = exact, mra = mra, lowrank = lowrank)



Linvs = list(truth = solve(t(chol(Sig0Mat))))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}



d = data.frame(truth = c(), exact = c(), mra = c())

for (iter in 1:maxiter) {

    cat(sprintf("%s iter %d\n", Sys.time(), iter))
    est = list()
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                     sig2 = sig2, smooth = smooth, range = range, locs = locs)

    L = Linvs[["truth"]]
    sse = 0
    for (t in 2:Tmax) {
        xs1 = XY$x[[t - 1]]
        xs2 = XY$x[[t]]
        w = xs2 - evolFun(xs1)
        sse = sse + t(w) %*% t(L) %*% L %*% w
    }
    est[["truth"]] = as.numeric(sse / ((N - 1) * (Tmax - 1)))
    

    for (name in names(approximations)) {

        cat(sprintf("%s \tsampling using %s\n", Sys.time(), name))
        app = approximations[[name]]
        smoothed = FFBS(app, XY$y, lik.params, prior_covparms, covparms, evolFun, Num_samples = 1)

        sse = 0
        Qinv = Linvs[[name]]
        states = smoothed$samples[[1]]
        for (t in 2:Tmax) {
            w = states[[t]] - evolFun(states[[t - 1]])
            sse = sse + (t(w) %*% Matrix::t(Qinv) %*% Qinv %*% w)
        }
        est[[name]] = as.numeric(sse/((N - 1)*(Tmax - 1)))
        
    }
        
    dloc = data.frame(truth = est[["truth"]], exact = est[["exact"]], mra = est[["mra"]], lr = est[["lowrank"]])
    d = rbind(d, dloc)
}

cat(sprintf("%s done!\n", Sys.time()))
cat(sprintf("Means of samples:\n"))
print(colMeans(d))

breaks = seq(range(d)[1], range(d)[2], length.out = 10)
ylim = c(0, 20)
oldpar = par(mfrow = c(4, 1))
hist(d$truth, xlab = "", ylim = ylim, xlim = range(d), breaks = breaks, main = "True data")
hist(d$exact, xlab = "", ylim = ylim, xlim = range(d), breaks = breaks, main = "Vecchia exact")
hist(d$mra, xlab = "", ylim = ylim, xlim = range(d), breaks = breaks, main = "Vecchia - hierarchical")
hist(d$lr, xlab = "sig2", ylim = ylim, xlim = range(d), breaks = breaks, main = "Low rank")
par(oldpar)
