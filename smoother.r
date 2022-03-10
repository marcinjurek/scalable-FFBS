## vecchia-smoothed means
smoothedMeans = function(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun) {

    preds = VecchiaFilter(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun,
                        saveUQ = c("L", "V"))
    n = length(Y[[1]])
    Tmax = length(Y)
    means = list()
    means[[Tmax]] = preds$preds[[Tmax]]
    
    for (t in (Tmax - 1):1) {

        E = evolFun(Matrix::Diagonal(n))
        L = preds$Ls[[t + 1]]
        V = preds$Vs[[t]]

        means[[t]] = means[[t + 1]] - preds$forecast[[t + 1]]
        means[[t]] = Matrix::solve(L, means[[t]], sparse = TRUE )
        means[[t]] = Matrix::solve( Matrix::t( L ), means[[t]], sparse = TRUE)
        means[[t]] = Matrix::t( E ) %*% means[[t]]
        means[[t]] = matrix(means[[t]], ncol = 1 )
        means[[t]] = Matrix::t(V) %*% means[[t]]
        means[[t]] = V %*% means[[t]]
        means[[t]] = preds$preds[[t]] + matrix(means[[t]], ncol = 1 )
        
    }
    return( list( means = means, preds = preds ) )
}
    


## exact Kalman Smoother
KalmanSmoother = function(Y, lik_params, prior_Cov, Q_cov,
                        evolFun, prior_mean = NULL, saveUQ = "L") {

    filtered = KalmanFilter(Y, lik_params, prior_Cov, Q_cov,
                            evolFun, prior_mean = prior_mean, saveUQ = saveUQ)
    n = length(Y[[1]])
    Tmax = length(Y)
    if (is.null(prior_mean)) {
        prior_mean = matrix(rep(0, n), ncol = 1)
    }
    E = evolFun(Matrix::Diagonal(n))
    
    smoothed = list()
    smoothedV = list()
    smoothed[[Tmax]] = filtered$preds[[Tmax]]
    smoothedV[[Tmax]] = solve(filtered$Vs[[Tmax]])
   
    for (t in (Tmax - 1):1) {
    
        L = filtered$Ls[[t + 1]]
        Linv = solve(L)
        V = filtered$Vs[[t]]
        mu_fi = filtered$preds[[t]]
        mu_fo = filtered$forecast[[t + 1]]
        J = V %*% Matrix::t(E) %*% Linv

        smoothed[[t]] = mu_fi + J %*% (smoothed[[t + 1]] - mu_fo)
        smoothedV[[t]] = V + J %*% (smoothedV[[t + 1]] - L) %*% Matrix::t(J)
     }
    
    return(list(smoothed = smoothed, smoothedV = smoothedV, filtered = filtered))
}
