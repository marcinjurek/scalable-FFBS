#!/usr/bin/env Rscript
cat(sprintf("%s Script started\n", Sys.time()))

source("getMatCov.r")
# we have to use the local version of the filter because we need
# the getMatCovFromFactor with the long WORD type and this currently
# doesn't work in GPvecchia
source("reusable-filter.r")
source("../gibbs-sampler/sig2-sample.r")
source("../gibbs-sampler/me-sample.r")
source("../advection.r")
#source("../aux-functions.r")
source("settings.r")
source("../simulate-data.r")
source("smoother.r")
source("sFFBS.r")


args = commandArgs(trailingOnly=TRUE)
init_covparms = c(SIG_02, RANGE, SMOOTH)
Qcovparms = c(SIG2, RANGE, SMOOTH)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, init_covparms)

## load data -----------------------------------
print(FILENAME)
TPW = readr::read_csv(FILENAME, col_types = cols()) %>%
    filter( lat > quantile(lat, 0.5), lon < quantile(lon, 0.5))

MEAN = TPW %>% summarize(mean(value, na.rm=TRUE)) %>% pull()
#MEAN = 0


locs = TPW %>% select(lon, lat) %>%
    data.matrix()
lat = sort(unique(locs[, "lat"]))
lon = sort(unique(locs[, "lon"]))[1:length(lat)]
locs = as_tibble(expand.grid(lon, lat))
colnames(locs) = c("lon", "lat")



Y = list()
Yf = list()
for (p in unique(TPW$day)) {
    if (p <= TMAX) {
        mu = TPW %>% dplyr::filter(day==p) %>% summarize(mean(value, na.rm=TRUE)) %>% pull()
        Yf[[p]] = TPW %>% dplyr::filter(day==p) %>% select(-day) %>% mutate(value = value - mu)
        Yf[[p]] = Yf[[p]] %>% right_join(locs, by = c("lon", "lat"))
        Yf[[p]] = Yf[[p]] %>%  arrange(lon, lat) %>% select(value) %>% pull()
        inds.obs = sample(1:length(Yf[[p]]), size = FRAC_OBS * length(Yf[[p]]))
        Y[[p]] = Yf[[p]]
        Y[[p]][-inds.obs] = NA
    }
}
nx = length(unique(locs$lon))
ny = length(unique(locs$lat))

save(list=c("Y", "Yf"), file="~/FFBS/data-application/results/TPW.Rdata")

name = args[1]
cat(sprintf("%s Setting up vecchia objects\n", Sys.time()))
if (name=="mra") {
    appr = GPvecchia::vecchia_specify(locs %>% data.matrix(), COND_SET_SIZE, conditioning = 'mra', mra.options = MRA_OPTIONS, verbose = TRUE)
} else if (name=="lrf") {
    appr = GPvecchia::vecchia_specify(locs %>% data.matrix(), COND_SET_SIZE, conditioning = 'firstm')
}



## predsLRF = VecchiaFilter(lrf, Y, lik.params, init_covparms, Qcovparms, evolFun,)
## predsLRF = predsMRA

cat(sprintf("%s Calculating the L matrices\n", Sys.time()))
Linvs = list()
md5 = digest::digest(appr, algo='md5')
filename = sprintf("%s/Lmat_%s.txt", WDIR, md5)
if (file.exists(filename)) {
    Linv = readMM(file=filename)
} else {
    matCov = GPvecchia::getMatCov(appr, covfun_d_0)
    L = GPvecchia::createL(appr, matCov)
    writeMM(L, file=filename)
    Linv = Matrix::solve(L, sparse = TRUE)
}




## filter ---------------------------------------


sigmas = rep(0, NSAMPLES)
mevars = rep(ME_VAR, NSAMPLES)
sigmas[1] = SIG2

lik_p = LIK_PARMS
    

for( sample_no in 2:NSAMPLES ) {

    cat(sprintf("%s Sampler: Working on sample %d\n", Sys.time(), sample_no))
    smoothed = sFFBS(appr, Y, lik_p, PRIOR_COVPARMS, covparams, evolFun, wDir = WDIR, Num_samples = 1, verbose=TRUE)
    states = smoothed$samples[[1]]
    mevars[sample_no] = mevar_sample(states, Y, A, B)
    sigmas[sample_no] = sig2_sample(states, evolFun, ALPHA, BETA, Linv)
    cat(sprintf("%s Sampler: sig2 = %f, me.var = %f\n", Sys.time(), sigmas[sample_no], mevars[sample_no]))
    covparams[1] = sigmas[sample_no]
    lik_p[["sigma"]] = sqrt(mevars[sample_no])
    save(list=c("states"), file=sprintf("~/FFBS/data-application/results/samples/sample-%s-%d.Rdata", name, sample_no))
    save(list=c("sigmas", "mevars"), file=sprintf("~/FFBS/data-application/results/results-%s.Rdata", name))

}
