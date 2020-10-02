FFBS = function(approx.name, obs, Nsamples = Nsamples, sig2draw = sig2){
    
    samples = list()
    for (sample.no in 1:Nsamples) {
    
        #cat(paste("Working on sample no", sample.no, "\n"))
    
        #cat("Simulating data\n")
        x0plus = matrix( RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1 )
        XYplus = simulate.xy(x0plus, evolFun, NULL, obs, lik.params, Tmax, sig2 = sig2draw, smooth = smooth, range = range, locs = locs)

        Y = mapply('-', obs, XYplus$y, SIMPLIFY = FALSE)

        browser()
        results = smoothedMeans(approx.name, Y, sig2draw)
       
        smeans = results[[ "means" ]]
        sample  = mapply('+', XYplus$x, smeans, SIMPLIFY = FALSE)

        samples[[sample.no]] = sample

    }
    return( list(filteringResults = results, samples = samples) )
}
