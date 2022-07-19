cat(sprintf("%s Script started\n", Sys.time()))
source("../aux-functions.r")
library(tidyverse)
#source("settings.r")

mevarsL = list()
sigmasL = list()
samples = list()

indobs = list()

load(file="~/FFBS/data-application/results/TPW.Rdata")  
TPWdata = readr::read_csv(file="/home/marcin/vecchiaFilter/data/katzfuss-dorit-data/MIRS.csv") %>%
        filter( lat > quantile(lat, 0.5), lon < quantile(lon, 0.5))

locs = TPWdata %>% select(lon, lat) %>%
    data.matrix()
lat = sort(unique(locs[, "lat"]))
lon = sort(unique(locs[, "lon"]))[1:length(lat)]
nx = length(unique(lon))
ny = length(unique(lat))
locs = expand.grid(lon, lat)


for (name in c("mra", "lrf")) {
    load(file=sprintf("~/FFBS/data-application/results/results-%s.Rdata", name))
    mevarsL[[name]] = mevars
    sigmasL[[name]] = sigmas
}
nsamples = length(mevarsL[["mra"]])
tmax = length(Y)
n = length(Y[[1]])



mu = rep(0, tmax)
for (p in unique(TPWdata$day)) {
    if (p <= tmax) {
        mu[[p]] = TPWdata %>% dplyr::filter(day==p) %>% summarize(mean(value, na.rm=TRUE)) %>% pull()
    }
}



M1 = list()
mulim = c(0, 0)
M2 = list()
sd = list()
sdlim = c(0, 0)
for (name in c("lrf", "mra")) {
    M1[[name]] = rep(list(rep(0, n)), tmax)
    M2[[name]] = rep(list(rep(0, n)), tmax)
    for (i in 2:nsamples) {
        load(file=sprintf("~/FFBS/data-application/results/samples/sample-%s-%s.Rdata", name, i))
        samples[[name]] = states
        for (t in 1:tmax) {
            M1[[name]][[t]] = M1[[name]][[t]] + as.numeric(states[[t]] + mu[t])
            M2[[name]][[t]] = M2[[name]][[t]] + (as.numeric(states[[t]] + mu[t])**2)
        }
    }
    sd[[name]] = list()
    for (t in 1:tmax) {
        M1[[name]][[t]] = M1[[name]][[t]] / (nsamples-1)
        mulim = range(c(mulim, M1[[name]][[t]]))
        M2[[name]][[t]] = M2[[name]][[t]] / (nsamples-1)
        sd[[name]][[t]] = sqrt(M2[[name]][[t]] - M1[[name]][[t]]**2)
        sdlim = range(c(sdlim, sd[[name]][[t]]))
    }
}



selectTs = c(2, 5, 8)

for (t in 1:tmax) {
    pdf(sprintf("~/FFBS/data-application/plots/TPW-%d.pdf", t), height=5, width=15)
    par(mfrow = c(1, 3))
    fields::quilt.plot(locs[, 1], locs[, 2], Y[[t]] + mu[t], nx=nx, ny=ny, zlim=mulim, main=sprintf("the data, t=%d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], M1[["mra"]][[t]], nx=nx, ny=ny, zlim=mulim, main=sprintf("scalable, t=%d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], M1[["lrf"]][[t]], nx=nx, ny=ny, zlim=mulim, main=sprintf("low-rank, t=%d", t))
    dev.off()    
}
    
#    M1[[name]] = list()
#    for (t in 1:tmax) {
#        M1[[t]] = sum(samples[[name]]
