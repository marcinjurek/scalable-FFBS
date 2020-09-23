source("advection.r")

N = 34**2
nx = ny = sqrt(N)
Tmax = 30

grid.oneside = seq(0,1,length = round(sqrt(N)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 

v = diffAdvVec2d(nx, ny = nx, height = 1, rnge = 4)
fields::quilt.plot( locs, as.numeric(v), nx = sqrt(N), ny = sqrt(N), main = "t=1" )

evolFun = function(X) evolAdvDiff(X, adv = 0.01, diff = 0.0)

for (t in 1:Tmax) {
  v = evolFun(v)
  if (t %% 3 == 0) {
    fields::quilt.plot( locs, as.numeric(v), nx = sqrt(N), zlim=c(-0.1, 1.1), ny = sqrt(N), main = paste("t=", t, sep = "") )  
  }
  
}