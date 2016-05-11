#' Locally Linear Embedding
#' @param data N x D matrix (N samples, D features)
#' @param dim integer; The target dimension.
#' @param k integer; Number of nearest neighbor.
#' @return Projected data; N x dim matrix
#' @details This is the exact same function as LLE in "RDRToolbox" except that the matrix inverse operation
#' is replaced by the pseudoinverse operation since the former often has numerical issues.
#' See "?RDRToolbox::LLE" for more details about usage.
#' @references Roweis, Sam T. and Saul, Lawrence K., "Nonlinear Dimensionality Reduction by Locally Linear Embedding", 2000
#' @examples
#' #Simulate data
#' sim_data <- swiss_roll(N = 600)
#' library(plotly)
#' p1 <- plotly_3D(sim_data); p1
#' lle_data <- LLE2(sim_data$data, dim = 2, k = 8)
#' p2 <- plotly_2D(lle_data, color = sim_data$colors); p2
#' @export
LLE2 <- function (data, dim = 2, k)
{
  if (missing(data))
    stop("data argument missing")
  else if (!is.matrix(data))
    stop("invalid argument: data argument is required to be a N x D matrix (N samples, D features)")
  if (!all(is.numeric(dim)) | dim < 1)
    stop("invalid argument: target dimension is required to be a positive integer value")
  if (missing(k))
    k = min(nrow(data), 5)
  else {
    if (k >= nrow(data))
      stop("invalid argument: more neighbours than samples")
    if (!is.numeric(k) | k <= 1)
      stop("invalid argument: neighbour parameter is required to be an integer value >= 2")
  }
  k = round(k)
  dim = round(dim)
  num_samples = nrow(data)
  num_features = ncol(data)
  message("Computing distance matrix ... ", appendLF = FALSE)
  d = t(rowSums(data^2))
  d = d[rep(1, each = num_samples), ]
  d = d + t(d) - 2 * data %*% t(data)
  diag(d) = 0
  d = sqrt(d)
  sort_idx = apply(d, 2, order)
  neighbours = sort_idx[2:(k + 1), ]
  message("done")
  message("Computing low dimensional emmbedding (using ", k,
          " nearest neighbours)... ", appendLF = FALSE)
  W = matrix(0, k, num_samples)
  for (i in 1:num_samples) {
    N = t(data[neighbours[, i], ]) - data[i, ]
    Cov = t(N) %*% N
    W[, i] = MASS::ginv(Cov) %*% rep(1, k)
    W[, i] = W[, i]/sum(W[, i])
  }
  M = diag(1, num_samples)
  for (i in 1:num_samples) {
    w = W[, i]
    n = neighbours[, i]
    M[i, n] = M[i, n] - t(w)
    M[n, i] = M[n, i] - w
    M[n, n] = M[n, n] + w %*% t(w)
  }
  eig_M = eigen(M)
  sweep(eig_M$vectors, 2, sqrt(colSums(eig_M$vectors^2)), "/")
  Y = eig_M$vectors[, (num_samples - 1):(num_samples - dim)] *
    sqrt(num_samples)
  message("done")
  return(Y)
}
