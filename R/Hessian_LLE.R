#' Hessian Local Linear Embedding
#' Written by David Donoho & Carrie Grimes, 2003.
#' Extracted from Todd Wittman's MANI: Manifold Learning Toolkit
#' Source: http://macs.citadel.edu/wittman/research.html
#' Port from Matlab to R by Jackson Kwok
HLLE <- function (X, k, d) {
  N = ncol(X)
  if (length(k) == 1) {
    kvec = rep(k, N)
  } else if (max(dim(k)) == N) {
    kvec = k
  }
  # Compute Nearest neighbors
  D1 = L2_distance(X, X, 1)
  dim = nrow(X)
  nind = matrix(numeric(nrow(D1) * ncol(D1)), nrow = nrow(D1))
  dp = d * (d + 1) / 2
  W = repmat(0, dp * N, N)

  if (mean(k) > d) {
    tol = 1e-3 # regularlizer in case constrained fits are ill conditioned
  } else {
    tol = 0
  }
  mse = numeric(N)
  for (i in 1:N) {
    tmp = D1[,i]
    or = order(tmp)
    #take k nearest neighbors
    nind[ or[2:(kvec[i] + 1)], i] = 1
    thisx = X[, or[2:(kvec[i] + 1)], drop = FALSE]
    #center using the mean
    mean_vec = apply(thisx, 1, mean)
    thisx = thisx - matrix(rep(mean_vec, kvec[i]), ncol = kvec[i])
    #compute local coordinates
    UDV = svd(thisx)
    Vpr = UDV$v
    V = Vpr[, 1:d]
    vals = UDV$d
    #Neighborhood diagnostics
    mse[i] = sum(tail(vals, -d))
    #build Hessian estimator
    Yi <- matrix(numeric(nrow(V) * d * (d+1) / 2), nrow = nrow(V), ncol = d * (d+1) / 2)
    ct = 0
    for (mm in 1:d) {
      startp = V[,mm]
      for (nn in seq_along(mm:d)) {
        indles = mm:d
        Yi[,ct+nn] = startp * V[,indles[nn]]
      }
      ct = ct + length(mm:d)
    }
    Yi = cbind(repmat(1, kvec[i], 1), V, Yi)
    #orthogonalize linear and quadratic forms
    mgs = modified_gram_schmidt(Yi)
    Yt = mgs$Q
    Orig = mgs$R

    Pii = t(Yt[ , (d+2):ncol(Yt)])
    # double check weights sum to 1
    for (j in 1:dp) {
      if (sum(Pii[j, ]) > 0.0001) {
        tpp = Pii[j, ] / sum(Pii[j, ])
      } else {
        tpp = Pii[j, ]
      }
      # fill weight matrix
      W[ (i - 1)*dp + j, or[2:(kvec[i] + 1)] ] = tpp
    }
  }
  #Compute eigenanalysis of W
  G = t(W) %*% W
  G = Matrix::Matrix(G, sparse = TRUE)
  G = as(G, "dgCMatrix")
  Y = t(RSpectra::eigs_sym(A = G, k = d + 1, sigma = -1e-14)$vectors[,1:d] * sqrt(N))

  #compute final coordinate alignment
  R = t(Y) %*% Y
  R2 = fnMatSqrtInverse(R)
  Y = Y %*% R2
  list(Y, mse)
}

#' Hessian Local Linear Embedding
#' @param X N x D matrix (N samples, D features).
#' @param k integer; Number of nearest neighbor.
#' @param d integer; The target dimension.
#' @return A list of two objects. The first is the projected data, the second is the mse.
#' @details Matlab codes were written by David Donoho & Carrie Grimes (2003) and
#' extracted from Todd Wittman's MANI: Manifold Learning Toolkit.
#' @references Donoho, D. L., & Grimes, C. (2003). Hessian eigenmaps: Locally linear embedding techniques for high-dimensional data. Proceedings of the National Academy of Sciences, 100(10), 5591-5596.
#' @references MANI: Manifold Learning Toolkit - http://macs.citadel.edu/wittman/research.html
#' @examples
#' #Simulate data
#' sim_data <- swiss_roll(N = 600)
#' library(plotly)
#' p1 <- plotly_3D(sim_data); p1
#' HLLE_data <- Hessian_LLE(sim_data$data, k = 8, d = 2)
#' p2 <- plotly_2D(HLLE_data$projection, color = sim_data$colors); p2
#' @export
Hessian_LLE <- function(X, k, d) {
  res = HLLE(t(X), k, d)
  list(projection = t(res[[1]]), mse = res[[2]])
}

#================================================================================
# utility function
#================================================================================
vector_norm <- function(v) {
  sqrt(sum(v^2))
}
# Modified Gram-Schmidt (Used by HLLE function.)
modified_gram_schmidt <- function(A){
  n <- ncol(A)
  R <- matrix(numeric(n^2), nrow = n)
  for (i in 1:n) {
    R[i,i] <- vector_norm(A[,i])
    A[,i] = A[,i] / R[i,i]
    if (i < n) {
      for (j in (i+1):n) {
        R[i,j] = sum(A[,i] * A[,j])
        A[,j] = A[,j] - R[i,j] * A[,i]
      }
    }
  }
  list(Q = A, R = R)
}
repmat <- function(v, nr, nc) {
  matrix(rep(v, nr * nc), nrow = nr)
}
# function to compute the inverse square root of a matrix
# reference: http://www4.stat.ncsu.edu/~li/software/SparseSDR.R
fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}
