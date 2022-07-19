rm(list=ls())
library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
source("data/process-data.r")


ndays = 6
m = 30 #5 10 20 40
frac_sample = 0.01
Ngrid = 40000

set.seed(1996)
params = read_csv("params_1_10.csv", col_types = cols())

xmin = 0; xmax = 0

data = list()
for( day.no in 1:ndays ){

    ## p = params %>% filter( day == day.no )
    ## TPW = read_day( day.no ) %>%
    ##     #filter( x < quantile(x, 0.1), y < quantile(y, 0.1) ) %>%
    ##     mutate( m = p$beta0 + p$beta1 * x )

    ## # training data
    ## data = TPW %>%
    ##     filter( !is.na(value) ) %>% sample_frac(frac_sample)
    ## locs = select( data, x, y ) %>% data.matrix()
    ## obs = select( data, value ) %>% pull()

    ## # test locations / data
    ## preds = TPW %>% filter( is.na(value) ) %>% sample_frac(frac_pred)
    ## cat("Remove TPW\n")

    ## locs.pred = select(preds, x, y) %>% data.matrix()
    
    ## cat("Calculating Vecchia approximation\n")    
    ## vecchia.approx = GPvecchia::vecchia_specify(locs,
    ##                                             conditioning = 'mra',
    ##                                             m = m,
    ##                                             locs.pred = locs.pred)

    ## cat("Calculating posterior\n")
    ## posterior = GPvecchia::calculate_posterior_VL(obs,
    ##                                               vecchia.approx,
    ##                                               likelihood_model = "gamma",
    ##                                               covparms = c(p$sig2, p$range, p$nu),
    ##                                               likparms = list("alpha" = p$a),
    ##                                               prior_mean = data$m,
    ##                                               return_all=TRUE)

    ## cat("Restructuring data for plotting\n")
    ## post = data %>% select(x, y) %>% add_column(value = posterior$mean, obs=TRUE )
    ## preds = mutate(preds, value = posterior$prediction$mu.pred + m) %>%
    ##     add_column( obs = FALSE ) %>%
    ##     bind_rows( post )

    cat( sprintf("data/latent/day_%d_ftrain_%.2f_Ngrid_%.2f_m_%d.csv\n", day.no, frac_sample, Ngrid, m) )
    data[[ day.no ]] = read_csv( sprintf("data/latent/day_%d_ftrain_%.2f_Ngrid_%.2f_m_%d.csv", day.no, frac_sample, Ngrid, m) ) %>%
        filter( row_number() > n() - Ngrid )
}


zlim = range(sapply(data, function(d) range(d$latent) ))


for( day.no in 1:ndays ){
    
    #if( day.no>1 ){

        #d = data[[ day.no ]] %>% mutate( w = latent - data[[ day.no - 1 ]]$latent ) %>% filter( is.na(value) )
        
        p = ggplot(data[[ day.no ]]) +
            geom_point( aes(x, y, fill = latent, stroke = 0 ), size = 2, shape = 21 ) +
        #p = ggplot(d) +
        #    geom_point( aes(x, y, fill = w, stroke = 0 ), size = 2, shape = 21 ) +
            scale_fill_gradientn(limits = zlim, colours = hcl.colors(10, palette = "Spectral", alpha = NULL, rev = FALSE, fixup = TRUE)) +
            ggtitle(sprintf("Latent field on Dec %d", day.no)) + xlab("longitude") + ylab("latitude")
        ggsave(sprintf("day_%d_ftrain_%.2f_Ngrid_%.2f_m_%d.png", day.no, frac_sample, Ngrid, m), p, width=7, height=7)
    #print(p)

    
            #p = ggplot(d) +
            #geom_point( aes(x, y, fill = value, stroke = 0 ), size = 2, shape = 21 ) +
        #p = ggplot(d) +
        #    geom_point( aes(x, y, fill = w, stroke = 0 ), size = 2, shape = 21 ) +
            #scale_fill_gradientn(limits = zlim, colours = hcl.colors(10, palette = "Spectral", alpha = NULL, rev = FALSE, fixup = TRUE)) +
            #ggtitle(sprintf("Latent field on Dec %d", day.no)) + xlab("longitude") + ylab("latitude")
        #ggsave(sprintf("obs_day_%d_ftrain_%.2f_Ngrid_%.2f_m_%d.png", day.no, frac_sample, Ngrid, m), p, width=7, height=7)
        #print(p)
    #}

}

