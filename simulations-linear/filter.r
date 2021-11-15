saveResults = function(preds.aux, L.tt, saveUQ){
    results = list(state = matrix(preds.aux$mean, ncol = 1))
    if("W" %in% saveUQ){
        results[["W"]] = preds.aux$W
    }
    if("L" %in% saveUQ) {
        results[["L"]] = L.tt
    }
    if("V" %in% saveUQ) {
        results[["V"]] = preds.aux$V
    }
    return( results )
}



forecastStep = function(appr, prior_mean, prior_covmodel, Qcovparms, evolFun) {
    cat(sprintf("\tForecast step\n"))
    n = nrow(appr$locs)
    covfun.d = function(D) GPvecchia::MaternFun(D, Qcovparms)
    cat(sprintf("\t\tCalculating Q matcov\n"))

    Qcovmodel = getMatCov(appr, covfun.d)
    E = evolFun(Matrix::Diagonal(n))

    mu.tt1 = as.numeric(E %*% prior_mean)
    Fmat = E %*% prior_covmodel
    cat(sprintf("\t\tCalculating forecast matcov\n"))
    priorMatCov = getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
    covmodel = Qcovmodel + priorMatCov
    Ltt1 = GPvecchia::createL( appr, covmodel )
    
    return(list(mean = mu.tt1, covmodel = covmodel, Ltt1 = Ltt1))
}



    
updateStep = function(obs, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ) {
    cat(sprintf("\tUpdate step\n"))
    loglik = NA    
    ## We need this because for some parameter values the posterior cannot be reasonably evaluated    
    tryCatch({
        cat(sprintf("\t\tCalculating posterior\n"))
        save(obs, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ, file="filter-test-data")
        preds.aux = GPvecchia::calculate_posterior_VL(obs, appr, prior_mean = forecast_mean,
                                                      likelihood_model = lik.params[["data.model"]],
                                                      covmodel = forecast_covmodel, covparms = NULL,
                                                      likparms = lik.params, return_all = TRUE,
                                                      max.iter=1000)
        if (!preds.aux$cnvgd) {
            warning("Posterior estimation did not converge")
        }
        L.tt = getLtt(appr, preds.aux)
        U = Matrix::t(L.tt)
        vars = as.numeric(sapply(split(U@x, cut(1:length(U@x), U@p, labels=FALSE)), function(v) sum(v**2)))
        mu.tt = matrix(preds.aux$mean, ncol = 1)
        preds.tt = saveResults(preds.aux, L.tt, saveUQ)
        msg = NULL
    },
    error = function(c) {
        msg = conditionMessage(c)
        logweight <<- -Inf
        loglik <<- -Inf
        cat(sprintf("%s\n", msg))
        stop(msg)
    })

    return(list(preds = preds.tt, vars = vars, msg = msg))
}




# filtering ---------------------

# We now extend this function to the case of multiple parameters.
# Here is a list of all parameters: likelihood params, c (multiplicative constant
# of the temporal evolution), covparms (constants in the matern model); that's the
# total of k=6.
KalmanFilter = function(appr, Y, lik.params, prior_covparms, prior_mean = NULL, saveUQ = "L"){

    if (is.null(appr)) {
        stop("No vecchia approximation provided")
    }

    n = nrow(appr$locsord)
    if (is.null(prior_mean)) {
        prior_mean = matrix(rep(0, n), ncol=1)
    }
    Tmax = length(Y)
    
    particles.all = list()
    preds.all = list()
    
    preds = list()
    Ls = list()
    Vs = list()
    forecast = list()

    # this can be modified to accomodate sampling particles for unknown
    # prior covariance models
    covfun.d = function(D) GPvecchia::MaternFun(D, prior_covparms)
    prior_covmodel = getL00(appr, covfun.d, locs)
    Qcovparms = c(sig2, range, smooth)

    for (t in 1:Tmax) {

        cat(sprintf("+++ Filtering for t = %d +++\n", t))
        forecasted = forecastStep(appr, prior_mean, prior_covmodel, Qcovparms, evolFun)
        updated = updateStep(Y[[t]], appr, forecasted$mean, forecasted$covmodel, lik.params, saveUQ)

        prior_mean = updated$preds$state
        prior_covmodel = updated$preds$L
        
        preds[[t]] = prior_mean#list()
        forecast[[t]] = forecasted$mean
        #preds[[t]][[1]] = list(state = prior_mean)
        Ls[[t]] = forecasted$Ltt1
        Vs[[t]] = updated$preds$V

    }
    
    return(list(preds = preds, forecast = forecast, Vs = Vs, Ls = Ls))
}



## KalmanFilter = function(approx, XY){

##     preds = list()

##     initCovfun.d = function(D) GPvecchia::MaternFun(D, prior_covparams)
##     prior_covmodel = getL00(approx, initCovfun.d)
    
##     for (t in 1:Tmax) {
      
##         cat(paste("filtering: t=", t, "\n", sep = ""))
##         yt = as.numeric(XY$y[[t]])
##         covmodel = GPvecchia::getMatCov(approx, covfun.d)
##         mu.tt1 = rep(0, n)
##         if(t > 1) {
##             Fmat = evolFun(L.tt)
##             M = GPvecchia::getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
##             covmodel = covmodel + M
##             mu.tt1 = evolFun( mu.tt )
##         }
##         #Ltt1 = GPvecchia::createL(approx, covparms, covmodel)      
        
##         preds.aux = GPvecchia::calculate_posterior_VL(yt, approx, prior_mean = mu.tt1,
##                                                       likelihood_model = data.model,
##                                                       covmodel = covmodel,
##                                                       likparms = lik.params,
##                                                       return_all = TRUE,
##                                                       covparms = covparms)
        
##         L.tt = getLtt(approx, preds.aux)
##         mu.tt = matrix(preds.aux$mean, ncol = 1)
##         preds[[t]] = list(state = mu.tt, var = L.tt)
##         #preds[[t]] = list(state = mu.tt, V = preds.aux$V, L = Ltt1, stateF = mu.tt1)
##     }    
##     return( preds )
## }
## filter = function(approx.name, Y){
  
##   approx = approximations[[ approx.name ]]
##   preds = list()

##   cat("Filtering\n")
    
##   covmodel = getMatCov(approx, covfun.d)
##   Ltt1 = GPvecchia::createL( approx, covparms, covmodel )
##   mu.tt1 = matrix( rep(0, n), ncol = 1 )
##   obs.aux = as.numeric(Y[[1]])

##   preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
##                                                 likelihood_model = data.model,
##                                                 covmodel = covmodel, covparms = covparms,
##                                                 likparms = lik.params, return_all = TRUE)
##   L.tt = getLtt(approx, preds.aux)
##   mu.tt = matrix(preds.aux$mean, ncol = 1)
##   preds[[1]] = list(state = mu.tt, V = preds.aux$V, L = Ltt1, stateF = mu.tt1)

  
##   if ( Tmax > 1 ) {
    
##     for (t in 2:Tmax) {
##       obs.aux = as.numeric(Y[[t]])
      
##       E = evolFun(Matrix::Diagonal(n))
##       Fmat = E %*% L.tt      
##       M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)        
##       M2 = getMatCov(approx, covfun.d)
##       covmodel = M1 + M2
##       Ltt1 = GPvecchia::createL( approx, covparms, covmodel )      
##       mu.tt1 = E %*% mu.tt

##       preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
##                                                       likelihood_model = data.model,
##                                                       covmodel = covmodel, covparms = covparms,
##                                                       likparms = lik.params, return_all = TRUE)

##       L.tt = getLtt(approx, preds.aux)
##       mu.tt = matrix(preds.aux$mean, ncol = 1)
##       preds[[t]] = list(state = mu.tt, V = preds.aux$V, L = Ltt1, stateF = mu.tt1)
##     }
##   }
##   return( preds )
## }
