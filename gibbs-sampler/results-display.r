writeParams = function(scores, resultsDir){

    avgRRMSPE = Reduce("+", lapply(scores, `[[`, 1))/length(scores)
    avgdLS = Reduce("+", lapply(scores, `[[`, 2))/length(scores)
    resultsAsString = list("===== avg. RRMSPE: ====\n")
    resultsAsString = c(resultsAsString, capture.output(print(avgRRMSPE)))
    resultsAsString = c(resultsAsString, "\n==== avg. dLS ====\n")
    resultsAsString = c(resultsAsString, capture.output(print(avgdLS)))
    resultsAsString = c(resultsAsString, "\n")
    resultsAsString = paste(resultsAsString, collapse = "\n")
    cat(resultsAsString)
    #output = paste(c(AllParamsAsString, resultsAsString), sep="\n")
    #writeLines(output, paste(resultsDir, "/logs/", hash, sep=""))

}



plot2Results = function(XY, means1, sds1, means2, sds2, path=".", names=NULL, diff=FALSE){

    Tmax = length( XY$x )
    range.truth = unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE)))
    range.data = unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))
    range.mean1 = unlist(lapply(means1, function(t) range(t, na.rm = TRUE)))
    range.mean2 = unlist(lapply(means2, function(t) range(t, na.rm = TRUE)))
    range.sd1 = unlist(lapply(sds1, function(t) range(t, na.rm = TRUE)))
    range.sd2 = unlist(lapply(sds2, function(t) range(t, na.rm = TRUE)))
    
    latent.obs.ratio = max(abs(range.truth))/max(abs(range.data))
    different.scales = latent.obs.ratio > 1.5 || latent.obs.ratio < 0.5
    
    range.latent = range(c(range.truth, range.mean1, range.mean2))
    if (!different.scales) {
        range.latent = range(c(range.latent, range.data))
        range.data = range.latent
    }

    range.sd = range( range.sd1, range.sd2 )
    
    for (t in 1:Tmax) {       
        if (t < 10) {
            number = paste("0", t, sep = "")  
        } else {
            number = t
        }
        
        truth = as.numeric(XY$x[[t]])
        mean1 = as.numeric(means1[[t]]) - diff*truth
        mean2 = as.numeric(means2[[t]]) - diff*truth
        sd1 = as.numeric(sds1[[t]])
        sd2 = as.numeric(sds2[[t]])
        
        png(paste(path, "/", number, ".png", sep = ""))#, width=8, height=8)
        defpar = par(mfrow = c(3, 2), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
        nna.obs = which(!is.na(XY$y[[t]]))
        fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim = range.data, nx = sqrt(n), ny = sqrt(n), main = "obs" )
        fields::quilt.plot( locs, truth, zlim = range.latent, nx = sqrt(n), ny = sqrt(n), main = "truth" )
        fields::quilt.plot( locs, sd1, zlim = range.sd, nx = sqrt(n), ny = sqrt(n), main = "sd 1" )
        fields::quilt.plot( locs, mean1, zlim = range.latent, nx = sqrt(n), ny = sqrt(n), main = "mean 1" )
        fields::quilt.plot( locs, sd2, zlim = range.sd, nx = sqrt(n), ny = sqrt(n), main = "sd 2" )
        fields::quilt.plot( locs, mean2, zlim = range.latent, nx = sqrt(n), ny = sqrt(n), main = "mean 2" )
        
        par(defpar)
        dev.off()
    }
}
 

   


plotResults = function(XY, means, sds, data.model){

    Tmax = length( XY$x )
    data1 = unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE)))
    data2 = unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))
    data3 = unlist(lapply(means, function(t) range(t, na.rm = TRUE)))
    
    if (data.model == 'poisson') {
        zrange = 1.1*range(c(data1, data3))
        obs.range = 1.1*range(data2)
    } else {
        obs.range = zrange = 1.1*range(c(data1, data2, data3))
    }
        
    for (t in 1:Tmax) {       
        if (t < 10) {
            number = paste("0", t, sep = "")  
        } else {
            number = t
        }
        png(paste(data.model, "/", number, ".png", sep = ""))#, width=8, height=8)
        defpar = par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
        nna.obs = which(!is.na(XY$y[[t]]))
        fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim = obs.range, nx = sqrt(n), ny = sqrt(n), main = "obs" )
        fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "truth" )
        fields::quilt.plot( locs, as.numeric(means[[t]]), zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "mean" )
        #fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "truth" )
        #fields::quilt.plot( locs, as.numeric(samples[[ sample.no ]][[t]]), zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "mean" )
        fields::quilt.plot( locs, as.numeric(sds[[t]]), nx = sqrt(n), ny = sqrt(n), main = "sd" )
        par(defpar)
        dev.off()
    }

}
