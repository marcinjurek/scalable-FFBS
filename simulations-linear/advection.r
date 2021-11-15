
######### simulate and plot the data #########
## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  # we assume that there is the same number of grid points
  # along each dimension
    if(is.null(dim(state))) {
        N = length(state)
    } else {
        N = dim(state)[1]
    }

    Ny = Nx = sqrt(N)
  
    dx = dy = 1/Nx
    D = diff/(dx**2)


    c1 = 1 + 2*(D + D) - adv*(1/dx + 1/dy)
    c2 = - D + adv*(1/dy)
    c3 = - D

    ## Let's consider the following layout on the grid
    ## . . . . .
    ## . . c . .
    ## . a x b .
    ## . . d . .
    ## . . . . .
    ## where the values "x" are given by d3, "d" - by d1,
    ## "a" - by d2, "b" - by d4 and "c" by d5.
    d = rep(c2, N-Nx)
    a = rep(c2, N-1)
    x = rep(c1, N)
    b = rep(c3, N-1)
    c = rep(c3,N-Nx)

    # Now notice that for k = 1:n, in each kn-th row there should
    # be no b and for every kn+1-st row there should be no a
    inds.b = Nx*(1:(Ny-1))
    inds.a = Nx*(1:(Ny-1))
    b[inds.b] = 0
    a[inds.a] = 0
    
    diags = list(d, a, x, b, c)
    E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)

    #browser()
    
    if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
        return( E %*% state )
    } else {
        return( as.numeric(E %*% as.matrix(state)) )
    }
}

  




diffAdvVec2d = function(nx, ny=nx, height=1, rnge=4){
  v = matrix(rep(0, nx*ny), ncol=ny)
  if((nx %% 2)==0) mid_x = nx/2 else mid_x = nx/2+1
  if((ny %% 2)==0) mid_y = ny/2 else mid_y = ny/2+1
  v[round(mid_x-rnge/2):round(mid_x+rnge/2), round(mid_y-rnge/2):round(mid_y+rnge/2)] = height
  
  return(matrix(as.numeric(v), ncol=1))
}


## ## define the temporal evolution function
## evolAdvDiff = function(state, adv=0, diff=0){
##   # we assume that there is the same number of grid points
##   # along each dimension
  
##   N = dim(state)[1]
##   Ny = Nx = sqrt(N)
  
##   dx = dy = 1/Nx
##   d = diff/(dx**2)
  
##   c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
##   c2 = - d + adv*(1/dy)
##   c3 = - d
  
##   diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
##   E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)
  
##   if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
##     return( E %*% state )
##   } else {
##     return( as.numeric(E %*% as.matrix(state)) )
##   }
## }



## diffAdvVec2d = function(nx, ny=nx, height=1, rnge=4){
##   v = matrix(rep(0, nx*ny), ncol=ny)
##   if((nx %% 2)==0) mid_x = nx/2 else mid_x = nx/2+1
##   if((ny %% 2)==0) mid_y = ny/2 else mid_y = ny/2+1
##   v[round(mid_x-rnge/2):round(mid_x+rnge/2), round(mid_y-rnge/2):round(mid_y+rnge/2)] = height
  
##   return(matrix(as.numeric(v), ncol=1))
## }
