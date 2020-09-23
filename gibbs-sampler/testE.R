n = 5**2
Nx = sqrt(n)
adv = 0.5





EMat = function(N, adv, diff){
  Ny = Nx = sqrt(N)

  dx = dy = 1/Nx
  d = diff/(dx**2)

  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = -d + adv*(1/dy)
  c3 = -d

  diags = list(rep(c2, N - Nx), rep(c2, N - 1), rep(c1, N), rep(c3, N - 1), rep(c3, N - Nx) )
  E = Matrix::bandSparse(N, k = c(-Nx, -1, 0, 1, Nx), diag = diags)
  return(E)
}


dx = 1/Nx
c1 = -1/dx - 1/dx
c2 = 1/dx
diags = list(rep(c2, n - Nx), rep(c2, n - 1), rep(c1, n))
A = Matrix::bandSparse(n, k = c(-Nx, -1, 0), diag = diags)

E = Matrix::Diagonal(n) + adv * A
