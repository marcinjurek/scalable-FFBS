sig2_sample = function(states, evolFun, alpha, beta, Qinv, numSamples=1) {

    n = length( states[[1]] )
    Tmax = length(states)
    new.alpha = alpha + n * (Tmax - 1) / 2
    new.beta = beta
    
    for (t in 2:Tmax) {
        w = states[[t]] - evolFun(states[[t - 1]])
        new.beta = new.beta + 0.5 * (t(w) %*% Matrix::t(Qinv) %*% Qinv %*% w)
    }

    draw = rinvgamma(numSamples, new.alpha, as.numeric(new.beta) )
    return( draw )
  
}
