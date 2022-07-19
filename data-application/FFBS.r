FFBS = function(approx, obs, lik.params, prior_covparams, Qcovparms,
                evolFun, wDir = NULL, Num_samples = 1, verbose = FALSE) {

    samples = list()
    n = nrow(approx$locsord)
    Tmax = length(obs)

    revord = order(approx$ord)
    locs = approx$locsord[revord, ]

    covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparms)
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L0 = GPvecchia::createL(appr, matCov)[approx$ord,]

    covfun_d_0 = function(D) GPvecchia::MaternFun(D, Qcovparms)
    matCovQ = GPvecchia::getMatCov(appr, covfun_d_0)
    LQ = GPvecchia::createL(appr, matCov)[approx$ord,]
    

    
    for (sample.no in 1:Num_samples) {

        if (verbose) {
            cat(sprintf("%s FFBS: Working on sample no %d\n", Sys.time(), sample.no))
            cat(sprintf("%s FFBS: Simulating data\n", Sys.time()))
        }

        #Sig0Model = RMwhittle(nu = prior_covparams[3], scale = prior_covparams[2], var = prior_covparams[1])
        #x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))
        
        
        XYplus = simulate.xy(x0, evolFun, NULL, obs, lik.params, Tmax, sig2 = Qcovparms[1],
                                  smooth = Qcovparms[3], range = Qcovparms[2], locs = locs)

        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        results = smoothedMeans(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun, wDir)

        smeans = results[[ "means" ]]
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)

        if (verbose) {
            cat(sprintf("%s FFBS: done working on sample %d\n", Sys.time(), sample.no))
        }
        samples[[sample.no]] = sample
        

    }
    return( list(filteringResults = results, samples = samples) )
}
