# scalable FFBS

This is the code and data used in obtaining simulation results for the paper "Scalable Spatio-Temporal Smoothing via Hierarchical Sparse Cholesky Decomposition". The data used in the application section can be found at:
https://github.com/katzfuss-group/vecchiaFilter/blob/publication/data-application/MIRS.csv

## How to reproduce the results
For a test run of any of the scripts we suggest first changing the settings to values that ensure quick computing time.

### Simulating the random field
Once the settings have been set in `simulations-linear/settings.r` run
```{r}
setwd("simulations-linear")
source("run-simulations.r")
```

### Calculating CRPS
CRP scores can be calculated using the `simulations-linear/calculate-crps.r` file. Once the results have been calculated the scores can be visualized using `plot-scores.r` script.


### Gibbs sampling
In order to run the Gibbs sampler, run
```{r}
setwd("gibbs-sampler")
source("run-simulations.r")
```


### Application to TPW
The application to TPW can be reproduced by running
```{r}
setwd("data-application")
source("run-sampler.r")
```
once the data file has been downloaded and the path set in the `data-application/settings.r` file.
