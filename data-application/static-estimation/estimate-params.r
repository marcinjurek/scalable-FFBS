library(Matrix)
suppressPackageStartupMessages(library(tidyverse))

set.seed(1996)
NDAYS = 9
M = 80
NOBS = 20000
SMOOTH = 1.5
FILENAME = "~/vecchiaFilter/data/katzfuss-dorit-data/MIRS.csv"
MAXITER = 500
PARMSFILE = "parms.csv"

params.names = c("m.e. var.", "sig2", "range", "mu")
params = matrix(rep(0, NDAYS*length(params.names)), ncol=length(params.names))
colnames(params) = params.names
params.file = "params.csv"

mus = rep(NA, NDAYS)

allTPW = readr::read_csv(FILENAME, col_types = cols()) #%>%
#    filter(lat > quantile(lat, 0.5), lon < quantile(lon, 0.5))

for (day_no in 1:NDAYS) {

    cat(sprintf("============== day %d =====================\n", day_no))

    TPW = allTPW %>% filter(day==day_no)
    
    mu = TPW %>% summarize(mean(value, na.rm=TRUE)) %>% pull()           
    mus[day_no] = mu

}
    
## cat(sprintf("Mean is: %f\n", mu))



for (day_no in 1:NDAYS) {

    cat(sprintf("============== day %d =====================\n", day_no))
    TPW = allTPW %>% filter(day==day_no) %>% filter(!is.na(value)) 
    n = nrow(TPW)

    if (n > NOBS) {
        TPW = TPW %>% sample_n(NOBS)
    }
    
    obs = TPW %>% pull(value)
    locs = TPW %>% select(lon, lat) %>% data.matrix()

    cat("Step 1, generating vecchia approximations\n")
    approx = GPvecchia::vecchia_specify(locs, conditioning = 'mra', m = M)

    ## Iterative method:  estimate a, then covparms, then a again
    cat("Step 2, optimizing parameters\n")

    initparms = c(1.5, 50, 1)
    cat("\tInitial parameter values:\n")
    print(initparms)

    t_start = Sys.time()

    hv_log_lik = function(x){
            covparms = c(exp(x[-1]), 1.5) ## sigma range smoothness
        lik_p = list(sigma=x[1], data.model = "gauss")
        #cat(sprintf("\t\tEvaluating covparms = (%.4f %.4f %.4f)\n", lik_p[["sigma"]], covparms[1], covparms[2]))
        
        ## Perform inference on latent mean with Vecchia Laplace approximation
        vll = GPvecchia::vecchia_laplace_likelihood(obs,
                                                    approx,
                                                    likelihood_model = "gauss",
                                                    covparms = covparms,
                                                    return_all = FALSE,
                                                    likparms = lik_p)
        #cat(sprintf("\t\tLikelihood for covparms = (%.4f %.4f %.4f): %.4f\n",
        #            lik_p[["sigma"]], covparms[1], covparms[2], vll))
        return(-vll)
    }


    x0 = log(initparms)
    vl = hv_log_lik(x0)

    opt_cont = list("trace" = 0,
                       "maxit" = MAXITER,
                       "reltol" = 1e-5,
                       "parscale" = x0)
    
    res = optim(x0, hv_log_lik, method = "Nelder-Mead", control = opt_cont)
    
    cat(sprintf("\tconvergence result: %s\n", res$convergence))
    cat(sprintf("\testimated parameter vector: %f, %f, %f\n", exp(res$par)[1], exp(res$par)[2], exp(res$par)[3]))

    t_end = Sys.time()
    time_dur = as.double(difftime(t_end, t_start, units = "mins"))
    cat(sprintf("Estimation took %.2f minutes\n", time_dur))

    params[day_no, ] = c(exp(res$par), mus[day_no])
    readr::write_csv(as.data.frame(params), file = PARMSFILE)
}
