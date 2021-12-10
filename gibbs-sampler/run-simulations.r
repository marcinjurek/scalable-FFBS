source('settings.r')

######### initialize #########
## generate grid of pred.locs
grid.oneside = seq(0, 1, length = round(sqrt(n)))
grid = expand.grid(grid.oneside, grid.oneside)
locs = as.matrix(grid) 


## set initial state
Sig0Model = RandomFields::RMwhittle(nu = smooth, scale = range, var = sig_02)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
Sig0Mat = RandomFields::RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)


RandomFields::RFoptions(storing = FALSE)
x0 = matrix( x0, ncol = 1)
zeros = matrix(rep(0, n), ncol = 1)


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm')
approximations = list(mra = mra)#, exact = exact)


covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparams)
Linvs = list(truth = solve(chol(Sig0Mat)))
for (name in names(approximations)) {
    appr = approximations[[name]]
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    Linvs[[name]] = solve(L, sparse = TRUE)
}
    


##### scalable FFBS
means = rep(0, max.iter)
for (iter in 1:max.iter) {
    print(iter)
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    XY_init = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

    all_draws = list(XY_init)
    all_sdraws = list()

    for (name in names(approximations)) {
        
        draws = rep(sig2, Nsamples)
        for (sample_no in 2:Nsamples) {
            
            sig2_draw = draws[sample_no - 1]
            covparams = c(sig2_draw, range, smooth)
            appr = approximations[[name]]
            sampled = FFBS(appr, XY$y, lik.params, prior_covparams, covparms, evolFun, Num_samples = 1, covparams = NULL)
            sampledStates = sampled$samples[[1]]

            if(Tmax>1) {
                sdraws = sig2.sample(sampledStates, evolFun, alpha, beta, Linvs[[name]], numSamples = 1000)
                draws[sample_no] = sdraws[1]
            }
            all_sdraws[[sample_no]] = sdraws
        }
        all_draws[[name]] = draws
    }
    means[iter] = mean(sdraws)
}

hist(means)
#plot(all_draws[["mra"]], pch = 16, col = "red")
#lines(all_draws[["exact"]], pch = 16, col = "black")

#par(mfrow = c(3, 3))
#for (sample_no in 2:(Nsamples-1)) {
#    if (sample_no %% 10 == 0) {
#        hist(all_sdraws[[sample_no]], main = sample_no)
#    }
#}

#title = paste("T=", Tmax, " , f=", frac.obs, ", me.var=", me.var, sep="")
#pdf("trajectory-sig2.pdf")
#plot(draws$sig2[20:Nsamples], type="l", xlab="sample no.", ylab="sig^2", main=title)
#dev.off()
