#' Laplacian Eigenmap.
#' Written by Belkin & Niyogi, 2002.
#' Extracted from Todd Wittman's MANI: Manifold Learning Toolkit
#' Source: http://macs.citadel.edu/wittman/research.html
#' Port from Matlab to R by Jackson Kwok
leigs <- function(data, type = 'nn', param, ne) {
  n = nrow(data)
  A = Matrix::Matrix(numeric(n*n), nrow = n)
  step = 100
  for (i1 in seq(1, n, step)) {
    i2 = i1 + step - 1
    if (i2 > n) {
      i2 = n
    }
    XX = data[i1:i2, ]
    dt = L2_distance(t(XX), t(data), 0)
    Z = t(apply(dt, 1, sort))
    I = t(apply(dt, 1, order))
    for (i in i1:i2) {
      ind = i - i1 + 1
      j = 2:(param + 1)
      A[i, I[ind, j]] = Z[ind, j]
      A[I[ind, j], i] = Z[ind, j]
    }
  }
  W = A
  W[W != 0] = 1
  D = apply(W, 1, sum)
  L = Matrix::Matrix(diag(D)) - W
  E = RSpectra::eigs_sym(L, ne, "SM") #default tolerence is 1e-10
  list(eigenvectors = E$vectors, eigenvalues = E$values)
}

#' Laplacian Eigenmap
#' @param X N x D matrix (N samples, D features).
#' @param k integer; Number of nearest neighbor.
#' @param d integer; The target dimension.
#' @return a list of two objects. The first is the projected data
#' (which are eigenvectors) and the second is the corresponding eigenvalues.
#' @details Matlab codes were written by Belkin & Niyogi (2001)
#' and extracted from Todd Wittman's MANI: Manifold Learning Toolkit.
#' @references Belkin, M., & Niyogi, P. (2001, December). Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering. In NIPS (Vol. 14, pp. 585-591).
#' @references MANI: Manifold Learning Toolkit - http://macs.citadel.edu/wittman/research.html
#' @examples
#' #Simulate data
#' sim_data <- swiss_roll(N = 600)
#' library(plotly)
#' p1 <- plotly_3D(sim_data); p1
#' LE_data <- Laplacian_Eigenmaps(sim_data$data, k = 8, d = 2)
#' p2 <- plotly_2D(LE_data$eigenvectors, color = sim_data$colors); p2
#' @export
Laplacian_Eigenmaps <- function(X, k, d) {
  proj = leigs(X, param = k, ne = d + 1)
  list(eigenvectors = proj$eigenvectors[, 1:d],
       eigenvalues = proj$eigenvalues[1:d])
}
