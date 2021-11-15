FFBS = function(approx, obs, lik.params, prior_covparams) {

    samples = list()
    n = nrow(approx$locsord)
    zeros = rep(0, n)
    for (sample.no in 1:Nsamples) {
    
        cat(paste("Working on sample no", sample.no, "\n"))
    
        cat("Simulating data\n")
        XYplus = simulate.xy(zeros, evolFun, NULL, obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
        w = XYplus$x[[2]] - XYplus$x[[1]]

        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        results = smoothedMeans(approx, Y, lik.params, prior_covparams)
    
        smeans = results[[ "means" ]]
        
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)
        #cat("done working on sample\n")
        samples[[sample.no]] = sample


    }
    return( list(filteringResults = results, samples = samples) )
}
