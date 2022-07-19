# initparms are (me.var, sig2, range)
fit_parms = function(initparms, approx, obs, maxiter) {


    hv_log_lik = function(x){
        covparms = c(exp(x[-1]), 1.5) ## sigma range smoothness
        lik_p = list(sigma=x[1], data.model = "gauss")
        cat(sprintf("\t\tEvaluating covparms = (%.4f %.4f %.4f)\n", lik_p[["sigma"]], covparms[1], covparms[2]))
        
        ## Perform inference on latent mean with Vecchia Laplace approximation
        vll = GPvecchia::vecchia_laplace_likelihood(obs,
                                                    approx,
                                                    likelihood_model = "gauss",
                                                    covparms = covparms,
                                                    return_all = FALSE,
                                                    likparms = lik_p)
        cat(sprintf("\t\tLikelihood for covparms = (%.4f %.4f %.4f): %.4f\n",
                    lik_p[["sigma"]], covparms[1], covparms[2], vll))
        return(-vll)
    }


    
    x0 = log(initparms)
    vl = vl_likelihood(x0)

    opt_cont = list("trace" = 0,
                       "maxit" = maxiter,
                       "reltol" = 1e-5,
                       "parscale" = x0)
    
    res = optim(x0, hv_log_lik, method = "Nelder-Mead", control = opt_cont)
    
    cat(sprintf("\tconvergence result: %s\n", res$convergence))
    cat(sprintf("\testimated parameter vector: %f, %f, %f\n", exp(res$par)[1], exp(res$par)[2], exp(res$par)[3]))


    return(exp(res$par))        
}

