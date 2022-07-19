library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
source("static-estimation.r")
source("../data/process-data.r")


ndays = 1
m = 50
n_sample = 10000


maxiter = 20
set.seed(1996)
niter.mean = 50



for( day.no in 1:ndays ){

    cat("====================================\n")
    TPW = read_day( day.no )
    n = nrow(TPW)
    
    #z = matrix( dplyr::select(TPW, value) %>% filter( !is.na(value) ) %>% pull(value), ncol=1 )
    #x = TPW %>% filter( !is.na(value) ) %>% pull(x)
    #y = TPW %>% filter( !is.na(value) ) %>% pull(y)
    z = TPW %>%
        dplyr::filter( x < -21.76952 ) %>%
        dplyr::filter( x < quantile(x, 0.2), y < quantile(y, 0.2) ) %>%
        #filter( !is.na(value) ) %>%
        select( value ) %>%
        pull()

    locs = TPW %>%
        dplyr::filter( x < -21.76952 ) %>%
        dplyr::filter( x < quantile(x, 0.2), y < quantile(y, 0.2) ) %>%
        #filter( !is.na(value) ) %>%
        select( x, y ) %>%
        data.matrix()

    
    # test on subset
    sub_idx = sample(length(z), n_sample, replace = FALSE)
    sub_locs = locs[sub_idx,]
    z_sub = z[sub_idx]

    
    #### trend estimation ####
    cat("Estimating mean\n")
    X = sub_locs[ !is.na(z_sub), ]
    zcomp = z_sub[ !is.na(z_sub) ]
    X[,1]= 1

    ## # Glm: IRLS method from taking deriv of llh
    beta = c(1, 0.001)
    for(i in c(1:niter.mean)){
      XB = X %*% beta
      W  = -Matrix::Diagonal(x = as.numeric(exp(-XB)*zcomp) )
      A  = exp(-XB)*zcomp-1
      U  = W %*% XB - A
      beta = solve( t(X) %*% W %*% X , t(X) %*% U)
    }
    #beta = c(7.546311, -0.083681)
    cat(sprintf("Mean coefficients: %f, %f\n", beta[1], beta[2]))
    XB = as.numeric(X %*% beta)




    cat("Generating vecchia approximations\n")
    vecchia.approx = GPvecchia::vecchia_specify(sub_locs, conditioning = 'mra', m = m)

    covparms = c(40, 5, 1.5)

    cat("Calculating likelihood\n")
    vll = GPvecchia::vecchia_laplace_likelihood(z_sub,
                                                vecchia.approx,
                                                likelihood_model="gamma",
                                                covparms = covparms,
                                                likparms = list(alpha=0.124979),
                                                prior_mean = as.numeric(XB))#,

    cat(sprintf("Likelihood is %f\n", vll))
}    
