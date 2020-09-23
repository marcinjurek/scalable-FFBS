a.sample = function(states, diffusion, mu.a, sig2.a, Qinv, sig2){
  
  n = length( states[[1]] )
  Nx = sqrt(n)
  num = mu.a/sig2.a
  denom = 1/sig2.a
  
  Qw.inv = Qinv/sig2
  c1 = 2; c2 = -1; c3 = -1
  diags = list(rep(c2, n - Nx), rep(c2, n - 1), rep(c1, n), rep(c3, n - 1), rep(c3,n - Nx) )
  D = Matrix::bandSparse(n, k = c(-Nx, -1, 0, 1, Nx), diag = diags)
  
  
  c1 = -1/sqrt(n) - 1/sqrt(n)
  c2 = 1/sqrt(n)
  diags = list(rep(c2, n - Nx), rep(c2, n - 1), rep(c1, n))
  A = Matrix::bandSparse(n, k = c(-Nx, -1, 0), diag = diags)
  
  for (t in 2:Tmax) {
    
    x.t = states[[ t ]] - (Matrix::Diagonal(n) + diffusion*D) %*% states[[ t - 1 ]]
    x.tt = A %*% states[[ t - 1 ]]
    
    num = num + as.numeric(t(x.tt) %*% Qw.inv %*% x.t)
    denom = denom + as.numeric(t(x.tt) %*% Qw.inv %*% x.tt)
  
  }
  
  new.mu = num/denom
  new.var = 1/denom

  cat(paste("new mu.a = ", new.mu, ", new sig2.a = ", new.var, "\n", sep = ""))
  return( abs(rnorm(1, new.mu, sd = sqrt(new.var))) )
  
}