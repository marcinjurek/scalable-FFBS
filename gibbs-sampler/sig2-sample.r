sig2.sample = function(states, evolFun, alpha, beta, Qinv) {
  
    n = length( states[[1]] )
    new.alpha = alpha + (Tmax - 1)*n/2
    
    new.beta = beta
    
    for (t in 2:Tmax) {
        w = states[[t]] - evolFun(states[[t - 1]])
        new.beta = new.beta + 0.5 * (t(w) %*% Qinv %*% w)
    }

    draw = rinvgamma(1, new.alpha, as.numeric(new.beta) )
    #est = t(w) %*% Qinv %*% w / (n-1)
    #ocat(paste("Sig2 =", draw, ", estimate = ", est, '\n'))
    return( draw )
  
}
