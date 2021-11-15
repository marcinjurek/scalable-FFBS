source('test-settings.r')

## generate grid of pred.locs
grid.oneside = seq(0, 1, length = round(sqrt(n)))
grid = expand.grid(grid.oneside, grid.oneside)
locs = as.matrix(grid) 

## set initial state
Sig0Model = RMwhittle(nu = smooth, scale = range, var = 1)
x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))

RFoptions(storing = FALSE)
Sig0Mat = RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)

## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = "firstm")
approximations = list(exact = exact, mra = mra)



PrecMat = solve(Sig0Mat)
varps = list(truth = PrecMat)

maxiter = 30
d = data.frame(truth = c(), exact = c(), mra = c())
for (iter in 1:maxiter) {
    cat(sprintf("======== iter %d ========\n", iter))
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                     sig2 = sig2, smooth = smooth, range = range, locs = locs)
    es = list()
    es[["truth"]] = XY$x[[2]] - XY$x[[1]]
    
    for (name in names(approximations)) {
        app = approximations[[name]]
        filtered = KalmanFilter(app, XY$y, lik.params, c(1, range, smooth))
        xf1 = filtered$preds[[1]]
        xf2 = filtered$preds[[2]]
        es[[name]] = xf2 - xf1
        if (name=="exact" & iter==1) {
            L00 = filtered$pvar
            varps[[name]] = solve(L00 %*% Matrix::t(L00))
        } else if(name=="mra") {
            L00 = filtered$pvar
            varps[[name]] = solve(L00 %*% Matrix::t(L00))
        }
               
    }

    est = list()
    for (name in names(es)) {
        e = es[[name]]
        est[[name]] = as.numeric(t(e) %*% varps[[name]] %*% e / (n - 1))
    }
    dloc = data.frame(truth = est[["truth"]], exact = est[["exact"]], mra = est[["mra"]])
    d = rbind(d, dloc)
}


## oldpar = par(mfrow = c(3, 1))
## zlim = c(0, 0)
## for (name in names(es)) {
##     e = es[[name]]
##     zlim = range(zlim, e)
##     fields::quilt.plot(locs[, 1], locs[, 2], e, nx = sqrt(n), ny = sqrt(n), main=name)
## }
## par(oldpar)
