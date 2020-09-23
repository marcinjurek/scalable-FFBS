sig2.sample = function(states, evolFun, alpha, beta, Qinv){
  
  n = length( states[[1]] )
  new.alpha = alpha + (Tmax - 1)*n/2
  
  new.beta = beta

  for (t in 2:Tmax) {
    v = states[[t]] - evolFun(states[[t - 1]])
    new.beta = new.beta + 0.5 * (t(v) %*% Qinv %*% v)
  }

  cat(paste("new.alpha = ", new.alpha, ", new.beta = ", new.beta, "\n", sep = ""))
  return( rinvgamma(1, new.alpha, as.numeric(new.beta) ) )
  
}