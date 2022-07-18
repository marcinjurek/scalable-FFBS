sig2_sample = function(states, evolFun, alpha, beta, Qinv, numSamples=1) {

    n = length( states[[1]] )
    Tmax = length(states)
    new.alpha = alpha + Tmax - 1
    new.betas = rep(beta, n)

    
    for (t in 2:Tmax) {
        w = states[[t]] - evolFun(states[[t - 1]])
        wtilde = Qinv %*% w
        new.betas = new.betas + 0.5 * t(wtilde) %*% wtilde
    }

    draws = rep(NA, n)
    for (j in 1:n) {
        draws[j] = rinvgamma(numSamples, new.alpha, new.betas[j])
    }
    return( draws )
  
}
