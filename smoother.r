revMat = function(mat) mat[nrow(mat):1,ncol(mat):1,drop = FALSE]

smoothedMeans = function(approx, Y, lik.params, prior_covparams, Qcovparms, evolFun) {

    preds = VecchiaFilter(appr, Y, lik.params, prior_covparms, Qcovparms, evolFun,
                        saveUQ = c("L", "V"))

    Tmax = length(Y)
    cat("Smoothing\n")
    means = list()
    means[[ Tmax ]] = preds$preds[[ Tmax ]]

    for (t in (Tmax - 1):1) {

        #cat(paste("Smoothing: t=", t, "\n", sep = ""))
        
        E = evolFun( Matrix::Diagonal( n ) )
        L = preds$Ls[[ t ]]
        V = preds$Vs[[ t ]]

        #means[[ t ]] = means[[ t + 1 ]] - preds$forecast[[ t + 1 ]]
        #means[[ t ]] = solve( L, means[[ t ]][ or ], sparse = TRUE )
        #means[[ t ]] = Matrix::t( E ) %*% solve( Matrix::t( L ), means[[ t ]], sparse = TRUE)
        #means[[ t ]] = matrix( rev( means[[ t ]] ), ncol = 1 )
        #means[[ t ]] = solve( Matrix::t( V ), solve( V, means[[ t ]], sparse = TRUE ), sparse = TRUE )
        #means[[ t ]] = preds$preds[[ t ]] + matrix( rev( means[[ t ]] )[ order( or ) ], ncol = 1 )

        means[[ t ]] = means[[ t + 1 ]] - preds$forecast[[ t + 1 ]]
        means[[ t ]] = solve(L, means[[ t ]], sparse = TRUE )
        means[[ t ]] = solve( Matrix::t( L ), means[[ t ]], sparse = TRUE)
        means[[ t ]] = Matrix::t( E ) %*% means[[ t ]]
        means[[ t ]] = matrix(means[[ t ]], ncol = 1 )
        means[[ t ]] = Matrix::t(V) %*% means[[ t ]]
        means[[ t ]] = V %*% means[[ t ]]
        means[[ t ]] = preds$preds[[ t ]] + matrix(means[[ t ]], ncol = 1 )
        
    }
    #cat("done smoothing\n")
    #return( means )

    return( list( means = means, preds = preds ) )
}
    


## KalmanSmoother = function(Y){

##     cat(paste("Filtering"))
##     filtered = KalmanFilter('exact', Y)
##     Tmax = length(Y)
##     or = approximations[['exact']]$ord
    
##     smoothed = list()
##     vInv = solve( filtered[[ Tmax ]]$V, sparse = TRUE )
##     sVar = revMat( t(vInv) %*% vInv )[ order(or), order(or) ]
##     smoothed[[ Tmax ]] = list(state = filtered[[ Tmax ]]$state, var = sVar)
    
##     E = evolFun( Matrix::Diagonal( n ) )[ or, or ]
  
##     for (t in (Tmax - 1):1) {
    
##         cat(paste("Smoothing: t=", t, "\n", sep = ""))
    
##         L = filtered[[ t ]]$L
##         L.ord = filtered[[ t ]]$L[ order( or ), order( or ) ]
##         V = filtered[[ t ]]$V
##         Linv = solve( L, sparse = TRUE )
##         Vinv = solve( V, sparse = TRUE )
        
##         C = revMat( Matrix::t( Vinv ) %*% Vinv %*% revMat( Matrix::t( E ) %*% Matrix::t( Linv ) %*% Linv ) )[ order( or ), order( or ) ]
    
##         smeans = smoothed[[ t + 1 ]]$state - filtered[[ t + 1 ]]$stateF        
##         smeans = Matrix::t(E) %*% Matrix::t(Linv) %*% Linv %*% smeans[ or ]
##         smeans = matrix( rev( smeans ), ncol = 1 )
##         smeans = Matrix::t(Vinv) %*% Vinv %*% smeans
##         smeans = filtered[[ t ]]$state + matrix( rev( smeans )[ order(or) ], ncol = 1 )
        
##         sVar = smoothed[[ t + 1 ]]$var - L.ord %*% Matrix::t(L.ord)
##         sVar = L.ord %*% Matrix::t(L.ord) + C %*% sVar %*% Matrix::t(C)
        
##         smoothed[[ t ]] = list( state = smeans, var = sVar )
##     }
  
##     return(list(smoothed = smoothed, filtered = filtered))
## }
