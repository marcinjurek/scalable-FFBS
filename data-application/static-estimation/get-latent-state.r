rm(list=ls())
library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
setwd("~/vecchiaFilter")
source("~/vecchiaFilter/aux-functions.r")
source("data-application/data/process-data.r")


ndays = 1
m = 30 #5 10 20 40
frac_sample = 1.0
#fraco_pred = 0.05

X_UPPER_LIMIT = -21.76952
LOC_QUANTILE = 0.1


set.seed(1996)
params = read_csv("~/vecchiaFilter/data-application/static-estimation/params_1_10.csv", col_types = cols())


grid = read_day( 1)  %>%
    dplyr::filter( x < X_UPPER_LIMIT ) %>%
    dplyr::filter( x > quantile(x, 1-LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) %>%
    data.matrix()
#    filter( x < quantile(x, q), y < quantile(y, q) ) %>%
#    summarise( range(x), range(y) ) %>%
#    data.matrix()
#x_seq = seq(from=limits[,'range(x)'][1], to=limits[,'range(x)'][2], length.out = Nx)
#y_seq = seq(from=limits[,'range(y)'][1], to=limits[,'range(y)'][2], length.out = Ny)
#grid = as.matrix(expand.grid(x = x_seq, y = y_seq))


day.no=1
for( day.no in 1:ndays ){

    p = params %>% filter( day == day.no )
    TPW = read_day( day.no ) %>%
        dplyr::filter( x < X_UPPER_LIMIT ) %>%
        dplyr::filter( x > quantile(x, 1-LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) %>%
        mutate(value = value/1000)
    
    # training data
    data = TPW %>% 
        filter( !is.na(value) ) %>% sample_frac(frac_sample)
    obs = select( data, value ) %>% pull()

    #preds = tibble(x = grid[,1], y = grid[,2], value = rep(NA, nrow(grid)), day = day.no)
    data = data %>% mutate( M = p$beta0 + p$beta1 * x )
    obs = select( data, value ) %>% pull()
    locs = select( data, x, y ) %>% data.matrix()

    # test locations / data
    #preds = TPW %>% filter( is.na(value) ) %>% sample_frac(frac_pred)
    #cat("Remove TPW\n")

    #locs.pred = grid#select(preds, x, y) %>% data.matrix()

    cat("Calculating Vecchia approximation\n")    
    vecchia.approx = GPvecchia::vecchia_specify(locs,
                                                m = m,
                                                conditioning = 'mra')#,
                                                #locs.pred = grid)

    cat("Calculating posterior\n")
    posterior = GPvecchia::calculate_posterior_VL(obs,
                                                  vecchia.approx,
                                                  likelihood_model = "gamma",
                                                  covparms = c(p$sig2, p$range, p$nu),
                                                  likparms = list("alpha" = p$a),
                                                  prior_mean = data$M,
                                                  return_all = TRUE)

    L.tt = getLtt(vecchia.approx, posterior)
    U = Matrix::t(L.tt)
    vars = as.numeric(sapply(split(U@x, cut(1:length(U@x), U@p, labels=FALSE)), function(v) sum(v**2)))

    
    obs_mean = exp(posterior$mean + vars/2)
    data = add_column( data, latent = obs_mean )
    #write_csv( data, sprintf("data/latent/day_%d_ftrain_%.2f_Ngrid_%.2f_m_%d.csv", day.no, frac_sample, nrow(grid), m) )
    p = ggplot(data) +
        geom_point( aes(x, y, fill = obs_mean, stroke = 0.001 ), size=1.0, shape=21 ) +
        scale_fill_gradientn(colours = hcl.colors(10, palette = "Spectral", alpha = NULL, rev = FALSE, fixup = TRUE)) +
        ggtitle(sprintf("Latent field on Dec %d", day.no)) + xlab("longitude") + ylab("latitude")
    ggsave(sprintf("day_%d_ftrain_%.2f_ftest_%.2f_m_%d.png", day.no, frac_sample, frac_sample, m), p)
    #print(p)
        
}
