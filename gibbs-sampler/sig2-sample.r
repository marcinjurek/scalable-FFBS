sig2.sample = function(states, evolFun, alpha, beta, Qinv, numSamples=1) {

    n = length( states[[1]] )
    new.alpha = n * (Tmax - 1) / 2 - 1
    new.beta = 0
    
    for (t in 2:Tmax) {
        w = states[[t]] - evolFun(states[[t - 1]])
        new.beta = new.beta + 0.5 * (t(w) %*% t(Qinv) %*% Qinv %*% w)
    }

    draw = rinvgamma(numSamples, new.alpha, as.numeric(new.beta) )
    return( draw )
  
}
