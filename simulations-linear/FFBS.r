FFBS = function(approx.name, obs){

    samples = list()
    for (sample.no in 1:Nsamples) {
    
        cat(paste("Working on sample no", sample.no, "\n"))
    
        XYplus = simulate.xy(zeros, evolFun, NULL, obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        results = smoothedMeans(approx.name, Y)
    
        smeans = results[[ "means" ]]
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)
        cat("done working on sample\n")
        samples[[sample.no]] = sample
    }
    smeans = getMeansFromSamples(samples)
    ssd = getSDFromSamples(samples)
    return( smoothedMeans = smeans, smoothedSD = ssd, filteringResults = results, samples = samples )
}