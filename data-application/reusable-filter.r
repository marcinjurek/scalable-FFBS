getLtt = function(vecchia.approx, preds){
  n = nrow(vecchia.approx$locsord)
  orig.order = order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse = TRUE)[seq(n, 1), ])[orig.order,]
  #browser()
  return(L.tt)
}


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



reusableForecastStep = function(appr, prior_mean, prior_covmodel, Qcovparms, evolFun, wDir = wDir) {

    n = nrow(appr$locs)
    E = evolFun(Matrix::Diagonal(n))
    mu.tt1 = as.numeric(E %*% prior_mean)

    hashlist = list(appr, Qcovparms, prior_covmodel)
    priorMD5 = digest::digest(hashlist, algo="md5")
    priorFile = sprintf("%s/%s", wDir, priorMD5)
    
    if (!file.exists(priorFile)) {
        covfun.d = function(D) GPvecchia::MaternFun(D, Qcovparms)
        Qcovmodel = getMatCov(appr, covfun.d)
        Fmat = E %*% prior_covmodel

        priorMatCov = getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
        covmodel = Qcovmodel + priorMatCov
        save(list=c("covmodel"), file=priorFile)
    } else {
        load(priorFile)
    }

    Ltt1 = GPvecchia::createL(appr, covmodel)
    
    return(list(mean = mu.tt1, covmodel = covmodel, Ltt1 = Ltt1))
}



    
reusableUpdateStep = function(obs, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ) {

    loglik = NA    
    ## We need this because for some parameter values the posterior cannot be reasonably evaluated    
    tryCatch({
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



VecchiaReusableFilter = function(appr, Y, lik.params, prior_covparms, Qcovparms, evolFun,
                                 prior_mean = NULL, saveUQ = c("L", "V"), wDir = ""){
    
    if (is.null(appr) | is(appr, "character")) {
        stop("No vecchia approximation provided")
    }

    n = nrow(appr$locsord)
    Tmax = length(Y)
    if (is.null(prior_mean)) {
        prior_mean = matrix(rep(0, n), ncol=1)
    }
    
    preds = list()
    Ls = list()
    Vs = list()
    fcasts = list()
    
    hashlist = list(appr, prior_covparms)
    priorMD5 = digest::digest(hashlist, algo="md5")
    priorFile = sprintf("%s/%s", wDir, priorMD5)

    if (!file.exists(priorFile)) {
        covfun.d = function(D) GPvecchia::MaternFun(D, prior_covparms)
        prior_covmat = getMatCov(appr, covfun.d)
        prior_covmodel = GPvecchia::createL(appr, prior_covmat)
        save(list=c("prior_covmodel"), file=priorFile)
    } else {
        covmodel = load(priorFile)
    }
        
    for (t in 1:Tmax) {

        cat(sprintf("%s Filter: Working on time %d\n", Sys.time(), t))
        fcast = reusableForecastStep(appr, prior_mean, prior_covmodel, Qcovparms, evolFun, wDir)
        updated = reusableUpdateStep(Y[[t]], appr, fcast$mean, fcast$covmodel, lik.params, saveUQ)

        prior_mean = updated$preds$state
        prior_covmodel = updated$preds$L

        preds[[t]] = prior_mean
        fcasts[[t]] = fcast$mean
        Ls[[t]] = fcast$Ltt1
        Vs[[t]] = updated$preds$L

    }
    
    return(list(preds = preds, forecast = fcasts, Vs = Vs, Ls = Ls))
}
