## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  # we assume that there is the same number of grid points
  # along each dimension
  
  N = dim(state)[1]
  Ny = Nx = sqrt(N)
  
  dx = dy = 1/Nx
  d = diff/(dx**2)
  
  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = -d + adv*(1/dy)
  c3 = -d
  
  diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
  E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)
  
  if (methods::is(state, 'matrix') || methods::is(state, 'sparseMatrix') || methods::is(state, 'Matrix')) {
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