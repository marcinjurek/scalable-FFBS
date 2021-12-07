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



forecastStep = function(appr, prior_mean, prior_covmodel, Qcovparms, evolFun) {

    n = nrow(appr$locs)
    covfun.d = function(D) GPvecchia::MaternFun(D, Qcovparms)

    Qcovmodel = GPvecchia::getMatCov(appr, covfun.d)
    E = evolFun(Matrix::Diagonal(n))

    mu.tt1 = as.numeric(E %*% prior_mean)
    Fmat = E %*% prior_covmodel

    priorMatCov = GPvecchia::getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
    covmodel = Qcovmodel + priorMatCov
    Ltt1 = GPvecchia::createL(appr, covmodel)
    
    return(list(mean = mu.tt1, covmodel = covmodel, Ltt1 = Ltt1))
}



    
updateStep = function(obs, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ) {

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



VecchiaFilter = function(appr, Y, lik.params, prior_covparms, Qcovparms, evolFun,
                        prior_mean = NULL, saveUQ = c("L", "V")){

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

    covfun.d = function(D) GPvecchia::MaternFun(D, prior_covparms)
    prior_covmat = GPvecchia::getMatCov(appr, covfun.d)
    prior_covmodel = GPvecchia::createL(appr, prior_covmat)

    for (t in 1:Tmax) {

        fcast = forecastStep(appr, prior_mean, prior_covmodel, Qcovparms, evolFun)
        updated = updateStep(Y[[t]], appr, fcast$mean, fcast$covmodel, lik.params, saveUQ)

        prior_mean = updated$preds$state
        prior_covmodel = updated$preds$L

        preds[[t]] = prior_mean
        fcasts[[t]] = fcast$mean
        Ls[[t]] = fcast$Ltt1
        Vs[[t]] = updated$preds$L
    }
    
    return(list(preds = preds, forecast = fcasts, Vs = Vs, Ls = Ls))
}

## exact Kalman Filter
KalmanFilter = function(Y, lik_params, prior_Cov, Q_cov,
                        evolFun, prior_mean = NULL, saveUQ = "L") {

    n = length(Y[[1]])
    Tmax = length(Y)
    if (is.null(prior_mean)) {
        prior_mean = matrix(rep(0, n), ncol=1)
    }
    
    preds = list()
    forecast = list()
    Ls = list()
    Vs = list()
    
    for (t in 1:Tmax) {
        ## forecast step
        obs = Y[[t]]
        E = evolFun(diag(n))
        f_mean = evolFun(prior_mean)
        f_Cov = E %*% prior_Cov %*% t(E) + Q_cov

        inds.obs = which(!is.na(obs))
        H = diag(n)[inds.obs, ]
        R = diag(length(inds.obs)) * (lik_params$sigma ** 2)
        K = f_Cov %*% t(H) %*% solve(H %*% f_Cov %*% t(H) + R)

        ##update step
        u_mean = f_mean + K %*% (obs[inds.obs] - H %*% f_mean)
        u_Cov = (diag(n) - K %*% H) %*% f_Cov

        prior_mean = u_mean
        prior_Cov = u_Cov
        
        preds[[t]] = u_mean
        forecast[[t]] = f_mean
        Ls[[t]] = f_Cov
        Vs[[t]] = u_Cov
    }
    return(list(preds = preds, forecast = forecast, Vs = Vs, Ls = Ls))
}
