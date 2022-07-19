cat(sprintf("%s Script started\n", Sys.time()))
source("../aux-functions.r")
library(tidyverse)
#source("settings.r")

mevarsL = list()
sigmasL = list()
samples = list()


load(file="~/FFBS/data-application/results/TPW.Rdata")  
TPWdata = readr::read_csv(file="/home/marcin/vecchiaFilter/data/katzfuss-dorit-data/MIRS.csv") %>%
        filter( lat > quantile(lat, 0.5), lon < quantile(lon, 0.5))

tmax=9
mu = rep(0, tmax)
for (p in unique(TPWdata$day)) {
    if (p <= tmax) {
        mu[[p]] = TPWdata %>% dplyr::filter(day==p) %>% summarize(mean(value, na.rm=TRUE)) %>% pull()
    }
}



indobs = list()
Yver = list()
load(file="~/FFBS/data-application/results/TPW.Rdata")
for (t in 1:length(Y)) {
    allobs = which(!is.na(Yf[[t]]))
    usedobs = which(!is.na(Y[[t]]))
    indobs[[t]] = setdiff(allobs, usedobs)
    Yver[[t]] = Yf[[t]][indobs[[t]]] + mu[t]
}
    




for (name in c("mra", "lrf")) {
    load(file=sprintf("~/FFBS/data-application/results/results-%s.Rdata", name))
    mevarsL[[name]] = mevars
    sigmasL[[name]] = sigmas
}
nsamples = length(mevarsL[["mra"]])


for (name in c("lrf", "mra")) {
    for (i in 2:nsamples) {
        load(file=sprintf("~/FFBS/data-application/results/samples/sample-%s-%s.Rdata", name, i))
        for (t in 1:length(Y)) {
            s = as.numeric(states[[t]])[indobs[[t]]]
            states[[t]] = s + mu[t]
        }       
        samples[[name]][[i - 1]] = states
        
    }
}


crps_mra = getSampleScores(Yver, samples[["mra"]])
crps_lrf = getSampleScores(Yver, samples[["lrf"]])
crps = crps_mra/crps_lrf


pdf(sprintf("~/FFBS/data-application/plots/CRPS-rates.pdf"), width=10, height=4)
ylim = range(crps)
plot(crps, col = "red", type = "l", ylim = ylim, ylab = "crps ratio", xlab = "time", lwd = 2)
dev.off()



pdf(sprintf("~/FFBS/data-application/plots/CRPS.pdf"), width=10, height=4)
ylim = range(c(crps_mra, crps_lrf))
plot(crps_mra, col="red", type="l", ylim=ylim, ylab="crps", xlab="time", lwd=2)
lines(crps_lrf, col = "blue", lwd = 2)
legend("topright", legend = c("scalable", "low-rank"), col = c("red", "blue"), lty = c(1, 1), lwd=2)
dev.off()
