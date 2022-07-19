## imports
library(Matrix)
suppressPackageStartupMessages(library(RandomFields))
library(invgamma)
suppressPackageStartupMessages(library(tidyverse))
`%>%` = dplyr::`%>%`


## settings -----------------------------------
COND_SET_SIZE = 55
#COND_SET_SIZE = 53
TMAX = 9
NSTEPS = 4
FRAC_OBS = 0.99
#FILENAME = "~/vecchiaFilter/data-application/data/katzfuss-dorit-data/MIRS.csv"
FILENAME = "~/vecchiaFilter/data/katzfuss-dorit-data/MIRS.csv"
#FILENAME = "~/FFBS/data-application/tests/simTPW.csv"
NSAMPLES = 100

MRA_OPTIONS = list(r = c(25, 10, 5, 5, 3, 2, 2), J = c(2, 2, 4, 4, 4, 4, 4, 4), M = 8)
#MRA_OPTIONS = list(r = c(25, 15, 15, 5, 5, 5, 3), J = 4, M = 7)
#MRA_OPTIONS = NULL

## evolution function -----------------------
C = 0.9
DIFFUSION = 0.000003
evolFun = function(X) {
    X = evolDiff(X, diff = DIFFUSION, nsteps = NSTEPS, rho = C)
    return(X)
}


## covariance parameters -------------------
#params = readr::read_csv("/home/marcin/vecchiaFilter/data-application/static-estimation/MIRS-Gaussian-params-smooth-1_5.csv", col_types = readr::cols())
params = readr::read_csv("~/vecchiaFilter/data-application/static-estimation/params_MIRS_1_5.csv", col_types = readr::cols())

SIG_02 = 83#
SIG2 = (1 - C) * SIG_02
RANGE = 1
SMOOTH = 1.5
PRIOR_COVPARMS = c(SIG_02, RANGE, SMOOTH)
covparams = c(SIG2, RANGE, SMOOTH)

## likelihood parameters -------------------
MEAN_COEFS = dplyr::select(params, mu) %>% data.matrix()
MEAN_COEFS[10,1] = sum(MEAN_COEFS[-10,1])/10

ME_VAR = 1.6
DATA_MODEL = "gauss"
LIK_PARMS = list(data.model = DATA_MODEL, sigma = sqrt(ME_VAR))

SIG2_MU = SIG2
SIG2_VAR = 0.000001 ** 2

ALPHA = (SIG2_MU ** 2)/(SIG2_VAR) + 2
BETA =  (SIG2_MU ** 3)/(SIG2_VAR) + SIG2_MU

ME_MU = ME_VAR
ME_V = 0.0000001 ** 2

A = (ME_MU ** 2)/(ME_V) + 2
B = (ME_MU ** 3)/(ME_V) + ME_MU

WDIR = "warehouse"


