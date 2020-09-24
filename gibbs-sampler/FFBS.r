FFBS = function(approx.name, obs, Nsamples = Nsamples, sig2draw = sig2){
    
    samples = list()
    for (sample.no in 1:Nsamples) {
    
        #cat(paste("Working on sample no", sample.no, "\n"))
    
        cat("Simulating data\n")
        XYplus = simulate.xy(zeros, evolFun, NULL, obs, lik.params, Tmax, sig2 = sig2draw, smooth = smooth, range = range, locs = locs)
        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)
        #print(sapply(Y, mean))
        results = smoothedMeans(approx.name, Y, sig2draw)
    
        smeans = results[[ "means" ]]
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)
        #cat("done working on sample\n")
        samples[[sample.no]] = sample
    }
    return( list(filteringResults = results, samples = samples) )
}