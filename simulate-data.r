## simulate y given x 
simulate.y = function(x, avail.obs, lik.params){

    if (is.numeric(x) & !is.matrix(x)) {
        x = matrix(x)
    }
    n = nrow(x)
    if( length( avail.obs ) == 1 ) {
        n.obs = round(n*avail.obs)
        obs.inds = sample(1:n, n.obs, replace = FALSE)
    } else {
        obs.inds = which( !is.na(avail.obs) )
        n.obs = length(obs.inds)
    }
  
    data.model = lik.params["data.model"]
    # simulate data
    if(data.model=='poisson'){
        y.obs = rpois(n.obs, exp(x[obs.inds]))
    } else if(data.model=='logistic'){
        y.obs = rbinom(n.obs,1,prob = exp(x[obs.inds])/(1+exp(x[obs.inds])))
    } else if(data.model=='gamma'){
        #default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)
        #z = rgamma(n.obs, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y[obs.inds]))
        y.obs = rgamma(n.obs, shape = lik.params[["alpha"]], rate = lik.params[["alpha"]]*exp(-x[obs.inds]))
        
    } else if(data.model == 'gauss'){
        y.obs = rnorm(n.obs, mean = x[obs.inds], sd=lik.params[["sigma"]])
               
    } else {
        print('Error: Distribution not implemented yet.')
    }
    y = rep(NA, n)
    y[obs.inds] = y.obs
    return(y)
}



## simulate x
simulate.xy = function(x0, E, Qc, avail.obs, lik.params, Tmax, seed=NULL, sig2=NULL, smooth=NULL, range=NULL, locs = NULL){

    if (!is.null(seed)) set.seed(seed)
    n = nrow(x0);
    x = list(); y = list()   

    #if (!is.null(Q) && any(Q)) {
    #    Qc = t(Matrix::chol(Q))
    #}
    for (t in 1:Tmax) {

        if (sig2>0 || !is.null(Qc) && sum(abs(Qc))>0) {
            if (!is.null(Qc)) {
                w =  Qc %*% matrix(rnorm(n), ncol = 1)
            } else {
                w = matrix(RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = smooth, scale = range, var = sig2),
                                                    x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1)
            }   
        } else {
            w = matrix(rep(0, n), ncol=1)
        }
        
        if (t==1) {
            x[[t]] = E(x0) + w
        } else {
            x[[t]] = E(x[[t - 1]]) + w
        }
        
        if( class( avail.obs ) == 'list' )  {
            y[[t]] = simulate.y(x[[t]], avail.obs[[ t ]], lik.params)
        } else {
            y[[t]] = simulate.y(x[[t]], avail.obs, lik.params)
        }
    } 

    return(list(x = x, y = y, x0 = x0))
}
