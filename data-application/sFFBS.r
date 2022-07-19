sFFBS = function(approx, obs, lik.params, prior_covparams, Qcovparms,
                evolFun, wDir = NULL, Num_samples = 1, verbose = FALSE) {

    samples = list()
    n = nrow(approx$locsord)
    Tmax = length(obs)

    revord = order(approx$ord)
    locs = approx$locsord[revord, ]

    if (verbose) {
        cat(sprintf("%s FFBS: Precomputing matrices\n", Sys.time()))
    }

    hash = digest::digest(list(approx, prior_covparams), algo="md5")
    filename = sprintf("%s/%s", WDIR, md5)
    if (file.exists(filename)) {
        L0 = Matrix::readMM(file=filename)
    } else {
        covfun_d_0 = function(D) GPvecchia::MaternFun(D, prior_covparams)
        matCov = GPvecchia::getMatCov(approx, covfun_d_0)
        L0 = GPvecchia::createL(approx, matCov)
        Matrix::writeMM(L0, file=filename)
    }

    
    hash = digest::digest(list(approx, Qcovparms), algo="md5")
    filename = sprintf("%s/%s", WDIR, md5)
    if(file.exists(filename)){
        LQ = Matrix::readMM(file=filename)
    } else {
        covfun_d = function(D) GPvecchia::MaternFun(D, Qcovparms)
        matCovQ = GPvecchia::getMatCov(approx, covfun_d)
        LQ = GPvecchia::createL(approx, matCov)
        Matrix::writeMM(LQ, file=filename)
    }

    
    for (sample.no in 1:Num_samples) {

        if (verbose) {
            cat(sprintf("%s FFBS: Working on sample no %d\n", Sys.time(), sample.no))
            cat(sprintf("%s FFBS: Simulating data\n", Sys.time()))
        }

        eps = matrix(rnorm(n), ncol=1)
        x0 = L0 %*% eps
        #x0 = matrix(RFsimulate(model = Sig0Model, x = locs[, 1], y = locs[, 2], spConform = FALSE))

        XYplus = simulate.xy(x0, evolFun, LQ, obs, lik.params, Tmax, locs = locs)

        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        results = smoothedMeans(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun, wDir, verbose=TRUE)

        smeans = results[[ "means" ]]
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)

        if (verbose) {
            cat(sprintf("%s FFBS: done working on sample %d\n", Sys.time(), sample.no))
        }
        samples[[sample.no]] = sample
        

    }
    return( list(filteringResults = results, samples = samples) )
}
