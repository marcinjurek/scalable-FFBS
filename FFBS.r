FFBS = function(approx, obs, lik.params, prior_covparams, Qcovparms,
                evolFun, Num_samples = 1, covparams = NULL) {

    if (!is.null(covparams)) {
        sig2 = covparams[1]
        range = covparams[2]
        smooth = covparams[3]
    }
    
    samples = list()
    n = nrow(approx$locsord)
    zeros = rep(0, n)
    for (sample.no in 1:Num_samples) {
    
        cat(paste("Working on sample no", sample.no, "\n"))
    
        cat("Simulating data\n")
        XYplus = simulate.xy(zeros, evolFun, NULL, obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        results = smoothedMeans(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun)
    
        smeans = results[[ "means" ]]
        
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)
        #cat("done working on sample\n")
        samples[[sample.no]] = sample


    }
    return( list(filteringResults = results, samples = samples) )
}
