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


Prec = solve(Sig0Mat)
XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax,
                 sig2 = sig2, smooth = smooth, range = range, locs = locs)

sse = 0
for (t in 2:Tmax) {
    w = XY$x[[t]] - evolFun(XY$x[[t - 1]])
    sse = sse + t(w) %*% Prec %*% w
}
mse = sse / ((N - 1) * (Tmax - 1))
cat(sprintf("%s marginal variance of the simulated process: %f\n", Sys.time(), as.numeric(mse)))
cat(sprintf("%s assumed marginal variance: %f\n", Sys.time(), sig2))

oldpar = par(mfrow = c(2, 3))
for (t in 1:Tmax) {
    v = XY$x[[t]]
    if (t %% 3 > 0) {
        next
    }
    fields::quilt.plot(locs, as.numeric(v), nx = sqrt(N), zlim = c(-4, 4),
                       ny = sqrt(N), main = paste("t=", t, sep = ""))
}
par(oldpar)
