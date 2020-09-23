########## smoothing ##########
revMat = function(mat) mat[nrow(mat):1,ncol(mat):1,drop = FALSE]

smoothedMeans = function(approx.name, Y){
   
    approx = approximations[[approx.name]]
    preds = filter(approx.name, Y)
    Tmax = length(Y)
    or = approx$ord
    
    cat("Smoothing\n")
    
    means = list()
    means[[ Tmax ]] = preds[[ Tmax ]]$state
    
    for (t in (Tmax - 1):1) {

        #cat(paste("Smoothing: t=", t, "\n", sep = ""))
        
        E = evolFun( Matrix::Diagonal( n ) )[ or, or ]
        L = preds[[ t ]]$L
        V = preds[[ t ]]$V
               
        means[[ t ]] = means[[ t + 1 ]] - preds[[ t + 1 ]]$stateF        
        means[[ t ]] = solve( L, means[[ t ]][ or ], sparse = TRUE )
        means[[ t ]] = Matrix::t( E ) %*% solve( Matrix::t( L ), means[[ t ]], sparse = TRUE)
        means[[ t ]] = matrix( rev( means[[ t ]] ), ncol = 1 )
        means[[ t ]] = solve( Matrix::t( V ), solve( V, means[[ t ]], sparse = TRUE ), sparse = TRUE )
        means[[ t ]] = preds[[ t ]]$state + matrix( rev( means[[ t ]] )[ order( or ) ], ncol = 1 )

    }

    #cat("done smoothing\n")
    #return( means )
    return( list( means = means, preds = preds ) )
}
    


KalmanSmoother = function(Y){

  cat(paste("Filtering"))
  filtered = filter('exact', Y)
  Tmax = length(Y)
  or = approximations[[ 'exact' ]]$ord
  
  smoothed = list()
  vInv = solve( filtered[[ Tmax ]]$V, sparse = TRUE )
  sVar = revMat( t(vInv) %*% vInv )[ order(or), order(or) ]
  smoothed[[ Tmax ]] = list(state = filtered[[ Tmax ]]$state, var = sVar)
  
  E = evolFun( Matrix::Diagonal( n ) )[ or, or ]
  
  for (t in (Tmax - 1):1) {
    
    cat(paste("Smoothing: t=", t, "\n", sep = ""))
    
    L = filtered[[ t ]]$L
    L.ord = filtered[[ t ]]$L[ order( or ), order( or ) ]
    V = filtered[[ t ]]$V
    Linv = solve( L, sparse = TRUE )
    Vinv = solve( V, sparse = TRUE )
    
    C = revMat( Matrix::t( Vinv ) %*% Vinv %*% revMat( Matrix::t( E ) %*% Matrix::t( Linv ) %*% Linv ) )[ order( or ), order( or ) ]
    
    smeans = smoothed[[ t + 1 ]]$state - filtered[[ t + 1 ]]$stateF        
    smeans = Matrix::t(E) %*% Matrix::t(Linv) %*% Linv %*% smeans[ or ]
    smeans = matrix( rev( smeans ), ncol = 1 )
    smeans = Matrix::t(Vinv) %*% Vinv %*% smeans
    smeans = filtered[[ t ]]$state + matrix( rev( smeans )[ order(or) ], ncol = 1 )
    
    sVar = smoothed[[ t + 1 ]]$var - L.ord %*% Matrix::t(L.ord)
    sVar = L.ord %*% Matrix::t(L.ord) + C %*% sVar %*% Matrix::t(C)
    
    smoothed[[ t ]] = list( state = smeans, var = sVar )
  }
  
  return(list(smoothed = smoothed, filtered = filtered))
}



########## filtering ##########
filter = function(approx.name, Y){
  
  approx = approximations[[ approx.name ]]
  preds = list()

  cat("Filtering\n")
  #cat(paste("filtering: t=1\n", sep = ""))
  #cat("\tCalculate covariance elements from function: ")
  #t0 = proc.time()
    
  covmodel = getMatCov(approx, covfun.d)
  Ltt1 = GPvecchia::createL( approx, covparms, covmodel )
   
  #t1 = proc.time()
  #cat(paste((t1 - t0)[3], "\n"))
  
  mu.tt1 = matrix( rep(0, n), ncol = 1 )
  obs.aux = as.numeric(Y[[1]])
  
  #cat("\tcalculate posterior: ")
  #t0 = proc.time()

  preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
                                                likelihood_model = data.model,
                                                covmodel = covmodel, covparms = covparms,
                                                likparms = lik.params, return_all = TRUE)

  #t1 = proc.time()
  #cat(paste((t1 - t0)[3], "\n"))
  #cat("\tstore results\n")
  #cat("\t\textract L.tt matrix\n")

  L.tt = getLtt(approx, preds.aux)

  #cat("\t\textract mean vector\n")

  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, V = preds.aux$V, L = Ltt1, stateF = mu.tt1)

  
  if ( Tmax > 1 ) {
    
    for (t in 2:Tmax) {
      
      #cat(paste("filtering: t=", t, "\n", sep = ""))

      obs.aux = as.numeric(Y[[t]])
      
      #cat("\tevolve the L.tt matrix:\n")
      #t0 = proc.time()
      #cat("\t\tbuild the evolution matrix\n")

      E = evolFun(Matrix::Diagonal(n))

      #cat("\t\tmultpily E by L\n")

      Fmat = E %*% L.tt
        
      #t1 = proc.time()
      #cat("\tCalculate covariance elements from factor: ")
      #t0 = proc.time()
      
      M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)

      #t1 = proc.time()
      #cat(paste((t1 - t0)[3], "\n"))
      #cat("\t... from function: ")
      #t0 = proc.time()
        
      M2 = getMatCov(approx, covfun.d)

      #t1 = proc.time()
      #cat(paste((t1 - t0)[3], "\n"))

      covmodel = M1 + M2
      Ltt1 = GPvecchia::createL( approx, covparms, covmodel )      
      mu.tt1 = E %*% mu.tt
        
      #cat("\tcalculate posterior: ")
      #t0 = proc.time()

      preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
                                                      likelihood_model = data.model,
                                                      covmodel = covmodel, covparms = covparms,
                                                      likparms = lik.params, return_all = TRUE)

      #t1 = proc.time()
      #cat(paste((t1 - t0)[3], "\n"))
      #cat("\tstore results\n")
      #cat("\t\textract L.tt matrix\n")

      L.tt = getLtt(approx, preds.aux)

      #cat("\t\textract mean vector\n")
      mu.tt = matrix(preds.aux$mean, ncol = 1)

      preds[[t]] = list(state = mu.tt, V = preds.aux$V, L = Ltt1, stateF = mu.tt1)
      
    }
    
  }
  return( preds )
}
