cat(sprintf("%s Script started\n", Sys.time()))

#source("getMatCov.r")
# we have to use the local version of the filter because we need
# the getMatCovFromFactor with the long WORD type and this currently
# doesn't work in GPvecchia
source("settings.r")

init_covparms = c(SIG_02, RANGE, SMOOTH)
Qcovparms = c(SIG2, RANGE, SMOOTH)
covfun_d_0 = function(D) GPvecchia::MaternFun(D, init_covparms)

## load data -----------------------------------
TPW = readr::read_csv(FILENAME, col_types = cols()) %>%
    filter( lat > quantile(lat, 0.0), lon < quantile(lon, 1.0))

locs = TPW %>% select(lon, lat) %>%
    filter( lat > quantile(lat, 0.0), lon < quantile(lon, 1.0)) %>%
    data.matrix()
lat = sort(unique(locs[, "lat"]))
lon = sort(unique(locs[, "lon"]))#[1:length(lat)]
locs = as_tibble(expand.grid(lon, lat))
colnames(locs) = c("lon", "lat")

#cat(sprintf("%s There are %d grid points\n". Sys.time(), nrow(locs)))
MEAN = TPW %>% summarize(mean(value, na.rm=TRUE)) %>% pull()

Y = list()
Yf = list()
means = rep(0, TMAX)
for (p in unique(TPW$day)) {
    if (p <= TMAX) {
        Yf[[p]] = TPW %>% dplyr::filter(day==p) %>% select(-day) %>% mutate(value = value - MEAN)
        Yf[[p]] = Yf[[p]] %>% right_join(locs, by = c("lon", "lat"))
        Yf[[p]] = Yf[[p]] %>%  arrange(lon, lat) %>% select(value) %>% pull()
        inds.obs = sample(1:length(Yf[[p]]), size = FRAC_OBS * length(Yf[[p]]))
        Y[[p]] = Yf[[p]]
        Y[[p]][-inds.obs] = NA
    }
}
nx = length(unique(locs$lon))
ny = length(unique(locs$lat))

x = select(locs, lon) %>% pull()
y = select(locs, lat) %>% pull()

par(mfrow = c(3, 4))
for (t in 1:TMAX) {
    fields::quilt.plot(x, y, Y[[t]], nx, ny, main = sprintf("t=%d", t), zlim = c(0, 40))
}

