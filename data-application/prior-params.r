### Pilot data
ranges = c(6.92, 6.898, 6.898, 5.873, 6.899)
as = c(0.1249, 0.1246, 0.1249, 0.1249, 0.1253)
sigs = c(27.12, 21.54, 21.54, 39.81, 21.54)
params = matrix(c(ranges, as, sigs), byrow=FALSE, ncol=3)
colnames(params) = c("range", "a", "sig2")


prior    = list(log.a = list(mean = mean(log(params[,"a"])), ## we assue that log(a) is normal
                             sd = diff(range(log(params[,"a"])))),
                
                ## c is normal
                c = list(mean = 0.8, sd = 0.3),
                
                ## we assume that inverse of sigma is gamma
                #inv.sig2 = list(shape = mean(1/params[,"sig2"])**2/var(1/params[,"sig2"]),
                #                rate = mean(1/params[,"sig2"])/var(1/params[,"sig2"])),
                ## .. or alternatively we use a log-normal distribution
                log.sig2 = list(mean = mean(log(params[,"sig2"])),
                                sd = diff(range(log(params[,"sig2"])))),

                
                ## we assume that log(range) is normal
                log.range = list(mean = mean(log(params[,"range"])),
                                 sd = diff(range(log(params[,"range"])))),
                
                ## we assume that nu is constant at 1.5
                log.nu = list(mean = log(1.5), sd = 0))


prop     = list(log.a = list( ## we assue that log(a) is normal
                    sd = 2 * diff(range(log(params[,"a"])))),
                
                ## c is normal
                c = list( sd = 0.6 ),
                
                ## we assume that inverse of sigma is gamma
                #inv.sig2 = list( shape = mean(1/params[,"sig2"])**2/var(1/params[,"sig2"]),
                #                 rate = mean(1/params[,"sig2"])/var(1/params[,"sig2"]) ),
                ## .. or alternatively we use a log-normal distribution
                log.sig2 = list(mean = mean(log(params[,"sig2"])),
                                sd = 2 * diff(range(log(params[,"sig2"])))),



                
                ## we assume that log(range) is normal
                log.range = list( sd = 2 * diff(range(log(params[,"range"])))),
                
                ## we assume that nu is constant at 1.5
                log.nu = list(sd = 0))



log.dist.eval = function( name, value, type ){

    if( type=="prior" ) {
        d = prior
    } else if( type=="proposal" ){
        d = prop
    } else {
        stop( "Wrong distribution type" )
    }
    
    if( name=="a" ) {
        dnorm( log(value), d$log.a$mean, d$log.a$sd, log=TRUE )
    } else if( name=="c" ) {
        dnorm( value, d$c$mean, d$c$sd, log=TRUE )
    } else if( name=="sig2" ) {
        dnorm( log(value), d$log.sig2$mean, d$log.sig2$sd, log=TRUE )
    } else if( name=="range" ) {
        dnorm( log(value), d$log.range$mean, d$log.range$sd, log=TRUE )
    } else if( name=="nu" ) {
        dnorm( log(value), d$log.nu$mean, d$log.nu$sd, log=TRUE )
    }
    
}


update.proposal = function( p ){

    prop$log.a$mean = log(p["a"])
    prop$c$mean = p["c"]
    prop$log.sig2$mean = log(p["sig2"])
    prop$log.range$mean = log(p["range"])
    prop$log.nu$mean = log(p["nu"])
    
}
