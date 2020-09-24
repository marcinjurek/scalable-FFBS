### This function is used mainly for testing.
### It takes the entire covariance matrix and creates
### a matrix of covariances
getL00 = function(vecchia.approx, covfun, locs){
  Sig.sel = GPvecchia::getMatCov(vecchia.approx, covfun(locs))
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = GPvecchia::ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j = inds, p = ptrs, x = vals, index1 = FALSE)
  ro = order(vecchia.approx$ord)
  return(Laux)
}



getLtt = function(vecchia.approx, preds){
  n = nrow(vecchia.approx$locsord)
  orig.order = order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse = TRUE)[seq(n, 1), ])[orig.order,]
  return(L.tt)
}

getCPRSscores = function(truth, smoothed){
  Tmax = length(truth)
  scores = rep(0, Tmax)
  samples.t = purrr:::transpose(samples)
  for (t in 1:Tmax) {
    mu = as.numeric(smoothed[[t]]$state)
    sig = as.matrix(smoothed[[t]]$var)
    scores[t] = scoringRules::crps_norm(as.numeric(truth[[t]]), mean = mu, sd = sig)
  }
  return(scores)
}

getScores = function(truth, samples = NULL, smoothed = NULL){
  if(!is.null(samples)){
    getSampleScores(truth, samples)
  } else if(!is.null(smoothed)){
    getCPRSscores(truth, smoothed)
  } else {
    stop("not enough data to evaluate the score")
  }
}

getSampleScores = function(truth, samples){
  Tmax = length(truth)
  scores = rep(0, Tmax)
  samples.t = purrr::transpose(samples)
  for (t in 1:Tmax) {
    dat = as.matrix(do.call(cbind, samples.t[[t]]))
    scores[t] = scoringRules::es_sample(as.numeric(truth[[t]]), dat)
  }
  return(scores)
}



getConfInt = function(preds, alpha){
  
  D = Matrix::diag(Matrix::solve(preds$W))
  mu = preds$state
  z = qnorm(1 - alpha/2)
  ub = mu + z*sqrt(D)
  lb = mu - z*sqrt(D)
  
  return(list(ub = ub, lb = lb))
}


getConfInt = function(preds, alpha){
  
  D = Matrix::diag(Matrix::solve(preds$W))
  mu = preds$state
  z = qnorm(1 - alpha/2)
  #browser()
  ub = mu + z*sqrt(D)
  lb = mu - z*sqrt(D)
  
  return(list(ub = ub, lb = lb))
}


## Calculate mean and standard deviation of samples at each grid point
getMeansFromSmoother = function(smoothed) {
  lapply(smoothed, function(t) t$state)
}
getSDFromSmoother = function(smoothed) {
  lapply(smoothed, function(t) diag(t$var))
}


getMeansFromFilter = function(filtered){
  lapply(filtered, function(t) t$state)
}

getSDFromFilter = function(filtered){
  
  getSD = function(dta){
    or = approximations[[1]]$ord
    V = dta$V
    Vinv = solve(V, sparse = TRUE)
    vars = rev(diag(Matrix::t(Vinv) %*% Vinv))[ order(or) ]
    sqrt(vars)
  }
  
  lapply(filtered, getSD)
}

getMeansFromSamples = function(samples){
  samples = purrr::transpose(samples)
  means = mapply( function(dta) Reduce('+', dta) / Nsamples, samples ) 
  return( means )
}

getSDFromSamples = function(samples){
  
  getSD = function(dta) {
    list.mean = Reduce('+', dta)/Nsamples
    list.squared.mean = Reduce("+", lapply( dta, "^", 2 )) / length( dta )
    list.variance = list.squared.mean - list.mean^2
    sqrt(list.variance)  
  }

    samples = purrr::transpose(samples)
  sds = mapply( getSD, samples )
  return( sds )
}
