cat(sprintf("%s Script started\n", Sys.time()))
source("settings.r")



params = list()
sigmas = list()

load(file="~/FFBS/data-application/results-mra.Rdata")
params[["mra"]] = params_field[["mra"]]
sigmas[["mra"]] = all_sigmas[["mra"]]
load(file="~/FFBS/data-application/results-lrf.Rdata")
params_field[["mra"]] = params[["mra"]]
all_sigmas[["mra"]] = sigmas[["mra"]]


init_covparms = c(SIG_02, RANGE, SMOOTH)
Qcovparms = c(SIG2, RANGE, SMOOTH)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, init_covparms)

## load data -----------------------------------
TPW = readr::read_csv(FILENAME, col_types = cols()) %>%
    filter( lat > quantile(lat, 0.5), lon < quantile(lon, 0.5))

locs = TPW %>% select(lon, lat) %>%
    data.matrix()
lat = sort(unique(locs[, "lat"]))
lon = sort(unique(locs[, "lon"]))[1:length(lat)]
locs = as_tibble(expand.grid(lon, lat))
colnames(locs) = c("lon", "lat")

if (!("lrf" %in% names(params_field))) {
    params_field[["lrf"]] = params_field[["mra"]]
    all_sigmas[["lrf"]] = all_sigmas[["mra"]]
}

if (!("mra" %in% names(params_field))) {
    params_field[["mra"]] = params_field[["lrf"]]
    all_sigmas[["mra"]] = all_sigmas[["lrf"]]
}


TMAX = length(params_field[["mra"]]$mu)
cat(sprintf("Tmax = %d\n", TMAX))

Y = list()
Yf = list()
for (p in unique(TPW$day)) {
    if (p <= TMAX) {
        Yf[[p]] = TPW %>% dplyr::filter(day==p) %>% select(-day)
        Yf[[p]] = Yf[[p]] %>% right_join(locs, by = c("lon", "lat"))
        Yf[[p]] = Yf[[p]] %>% arrange(lon, lat) %>% select(value) %>% pull()
        inds.obs = sample(1:length(Yf[[p]]), size = FRAC_OBS * length(Yf[[p]]))
        Y[[p]] = Yf[[p]]
        Y[[p]][-inds.obs] = NA
    }
}

nx = length(unique(locs$lon))
ny = length(unique(locs$lat))




MRF = LRF = sdMRF = sdLRF = list()
for (t in 1:TMAX) {
    MRF[[t]] = as.numeric(params_field[["mra"]][["mu"]][[t]])
    LRF[[t]] = as.numeric(params_field[["lrf"]][["mu"]][[t]])
    sdMRF[[t]] = sqrt(as.numeric(params_field[["mra"]][["sig2"]][[t]]))
    sdLRF[[t]] = sqrt(as.numeric(params_field[["lrf"]][["sig2"]][[t]]))
}




## plot results ------------------------
sdlim = c(1, 1)
xlim = c(1)
for (name in c("mra", "lrf")) {
    for (t in 1:TMAX) {
        xlim = range(c(xlim, Yf[[t]], Y[[t]], MRF[[t]], LRF[[t]]), na.rm=TRUE)
        sdlim = range(c(sdlim, sdMRF[[t]], sdLRF[[t]]), na.rm=TRUE)
    }
}

for (t in 1:TMAX) {
    #pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=5, height=20)
    pdf(sprintf("~/FFBS/data-application/plots/fields-%d.pdf", t), width=15, height=10)
    oldpar = par(mfcol = c(2, 3))
    
    fields::quilt.plot(locs$lon, locs$lat, Yf[[t]], zlim = xlim, nx = nx, ny = ny, main = sprintf("full data, t=%d", t))
    fields::quilt.plot(locs$lon, locs$lat, MRF[[t]], zlim = xlim, nx = nx, ny = ny, main = "predictions HV")
    fields::quilt.plot(locs$lon, locs$lat, sdMRF[[t]], zlim = sdlim, nx = nx, ny = ny, main = "sd HV")
    fields::quilt.plot(locs$lon, locs$lat, Y[[t]], zlim = xlim, nx = nx, ny = ny, main = sprintf("observations, t=%d", t))
    fields::quilt.plot(locs$lon, locs$lat, LRF[[t]], zlim = xlim, nx = nx, ny = ny, main = "predictions LRF")
    fields::quilt.plot(locs$lon, locs$lat, sdLRF[[t]], zlim = sdlim, nx = nx, ny = ny, main = "sd LRF")

     par(oldpar)
     dev.off()
}


pdf(sprintf("~/FFBS/data-application/plots/sigma-samples.pdf"), width=15, height=10)
siglim = range(as.numeric(unlist(all_sigmas)))
plot(all_sigmas[["mra"]], ylim = siglim, col="red", type = "l")
lines(all_sigmas[["lrf"]], col="blue", type = "l")
dev.off()


Ytr = list()
for (t in 1:TMAX) {
    ind.obs = which(!is.na(Y[[t]]))
    ind.pred = setdiff(which(!is.na(Yf[[t]])), ind.obs)
    Ytr[[t]] = Yf[[t]][ind.pred]
}

local_scores[[name]] = getScores(Ytr, smoothedSamples)



pdf(sprintf("~/FFBS/data-application/plots/MSE-rates.pdf"), width=15, height=10)
ylim = range(c(MSE_HV, MSE_LR))
plot(MSE_HV[1:9], col = "red", type = "l", ylim = ylim, ylab = "RMSE", xlab = "time", lwd = 2)
lines(MSE_LR[1:9], col = "blue", lwd = 2)
legend("topleft", legend = c("HV", "LR"), col = c("red", "blue"), lty = c(1, 1))
dev.off()
