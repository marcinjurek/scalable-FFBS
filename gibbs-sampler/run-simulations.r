source('/home/marcin/MRS/gibbs-sampler/settings.r')

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
Sig0Model = RandomFields::RMmatern(nu = smooth, scale = range)
x0 = RandomFields::RFsimulate(model = Sig0Model, x = locs[,1], y = locs[,2], spConform = FALSE)
Sig0Mat = RandomFields::RFcovmatrix(Sig0Model, distances = dist(locs), dim = 2)


RandomFields::RFoptions(storing = FALSE)
x0 = matrix( x0, ncol = 1)
zeros = matrix(rep(0, n), ncol = 1)


## define Vecchia approximation
#mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
mra = GPvecchia::vecchia_specify(locs, n-1, conditioning = 'firstm')
approximations = list(mra = mra)#, exact = exact)

Qinv = solve(Sig0Mat)

##### scalable FFBS
for (iter in 1:max.iter) {
    cat(paste("Starting iteration ", iter, "\n", sep = ""))

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    ###### naive Kalman Smoother
    #KSsmooth = KalmanSmoother(XY$y)

    sig2draw = sig2#runif(1)#rinvgamma(1, alpha, beta)
    cat(paste("Sig2 =", sig2draw, '\n'))
        
    #draws = list(state = list(), sig2 = rep(0, Nsamples))
        
    for (sample.no in 1:Nsamples) {
        cat(paste("Sample no", sample.no, "\n"))
        smoothingResults = FFBS('mra', XY$y, Nsamples = 1, sig2draw)
        sampledStates = smoothingResults$samples[[1]]
        #draws$state[[ sample.no ]] = sampledStates
        if(Tmax>1){
            sig2draw = sig2.sample(sampledStates, evolFun, alpha, beta, Qinv)
            #draws$sig2[ sample.no ] = sig2draw
        }
        v = chol(Qinv) %*% sampledStates[[1]]
        print(paste('x0 var:', t(v)%*%v/(n-1)))
    }
}


#title = paste("T=", Tmax, " , f=", frac.obs, ", me.var=", me.var, sep="")
#pdf("trajectory-sig2.pdf")
#plot(draws$sig2[20:Nsamples], type="l", xlab="sample no.", ylab="sig^2", main=title)
#dev.off()
