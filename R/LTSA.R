#' Local Tangent Space Alignment (LTSA)
#' Written by Zhenyue Zhang & Hongyuan Zha, 2004.
#' Extracted from Todd Wittman's MANI: Manifold Learning Toolkit
#' Source: http://macs.citadel.edu/wittman/research.html
#' Port from Matlab to R by Jackson Kwok
LTSA <- function(data, d, K, NI) {
  m <- nrow(data)  #number of features
  N <- ncol(data)
  # Step 0:  Neighborhood Index
  if (missing(NI)) {
    if (length(K) == 1) {
      K = rep(K, N)
    }
    NI = array(list(), N)
    if (m > N) {
      a = apply(data^2, 2, sum)
      dist2 = sqrt( matrix(rep(a, N), ncol = N) +
                    matrix(rep(a, N), nrow = N, byrow = T) -
                    2 * (t(data) %*% data) )
      for (i in 1:N) {
        # Determine ki nearest neighbors of x_j
        J = order(dist2[,i])
        Ii = J[1:K[i]]
        NI[[i]] = Ii
      }
    } else {
      for (i in 1:N) {
        # Determine ki nearest neighbors of x_j
        x = data[ ,i]
        ki = K[i]
        dist2 = apply( (data - matrix(rep(x, N), ncol = N))^2, 2, sum )
        J = order(dist2)
        Ii = J[1:ki]
        NI[[i]] = Ii
      }
    }
  } else {
    K = numeric(N)
    for (i in 1:N) {
      K[i] = length(NI[[i]])
    }
  }

  #Step 1:  local information
  BI = array(list(), N)
  for (i in 1:N) {
    # Compute the d largest right singular eigenvectors of the centered matrix
    Ii = NI[[i]]
    ki = K[i]

    tmp = apply(data[, Ii], 1, mean)
    Xi = data[,Ii] - matrix(rep(tmp, ki), ncol = ki)
    W = t(Xi) %*% Xi
    W = (W + t(W)) / 2

    tmp = Matrix::Schur(W)
    Vi = tmp$Q
    Si = tmp$T
    Ji = order(-diag(Si))
    Vi = Vi[,Ji[1:d]]

    # construct Gi
    Gi = cbind(matrix(rep(1 / sqrt(ki), ki), nrow = ki, byrow = TRUE), Vi)
    # compute the local orthogonal projection Bi = I-Gi*Gi'
    # that has the null space span([e,Theta_i^T]).
    BI[[i]] = diag(rep(1, ki)) - Gi %*% t(Gi)
  }
  B = Matrix::Matrix(diag(rep(1, N)))
  for (i in 1:N) {
    Ii = NI[[i]]
    B[Ii,Ii] = B[Ii,Ii] + BI[[i]]
    B[i,i] = B[i,i] - 1
  }
  B = (B + t(as.matrix(B))) / 2

  tmp = RSpectra::eigs_sym(B, d+2, sigma = -1e-14)
  U = tmp$vectors
  D = tmp$values
  lambda = D
  J = order(abs(lambda))
  U = U[ ,J]
  lambda = lambda[J]
  T = U[, 2:(d+1)]
  T
}

#' Local Tangent Space Alignment
#' @param X N x D matrix (N samples, D features).
#' @param k integer; Number of nearest neighbor.
#' @param d integer; The target dimension.
#' @return N x d matrix; the projected data.
#' @details Matlab codes were written by Zhenyue Zhang & Hongyuan Zha (2004)
#' and extracted from Todd Wittman's MANI: Manifold Learning Toolkit.
#' @references Zhang, Z., & Zha, H. (2004). Principal Manifolds and Nonlinear Dimensionality Reduction via Tangent Space Alignment. SIAM Journal on Scientific Computing, 26(1), 313-338.
#' @references MANI: Manifold Learning Toolkit - http://macs.citadel.edu/wittman/research.html
#' @examples
#' #Simulate data
#' sim_data <- swiss_roll(N = 600)
#' library(plotly)
#' p1 <- plotly_3D(sim_data); p1
#' LTSA_data <- Local_TSA(sim_data$data, k = 8, d = 2)
#' p2 <- plotly_2D(LTSA_data, color = sim_data$colors); p2
#' @export
Local_TSA <- function(X, k, d) {
  LTSA(t(X), d, k)
}
