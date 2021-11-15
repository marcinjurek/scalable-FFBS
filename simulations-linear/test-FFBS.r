source('settings.r')

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 

## set initial state
Sig0Model = RandomFields::RMwhittle(nu = smooth, scale = range, var = 1)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
x0 = matrix(x0)
RandomFields::RFoptions(storing = FALSE)

## function to calculate the variance estimate
Sig0Mat = RandomFields::RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)
invCorr = solve(Sig0Mat)
est_var = function(x) (t(x) %*% invCorr %*% x) / (n - 1)

## vecchia approximations
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm')
approximations = list(mra = mra, exact = exact)

## simulation
max.iter = 1
d = rep(0, max.iter)
for (i in 1:max.iter) {
    cat(sprintf("iteration %d\n", i))
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    predsE = filter('exact', XY$y)
    predsMRA = filter('mra', XY$y)
    eE = est_var(predsE[[2]]$state - predsE[[1]]$state)[1, 1]
    eMRA = est_var(predsMRA[[2]]$state - predsMRA[[1]]$state)[1, 1]
    d[i] = eMRA - eE
}

hist(d)



