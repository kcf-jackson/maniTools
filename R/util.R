#' Scale data to zero mean and unit variance
#' @param my_data matrix; rows correspond to records, columns correspond to features.
#' @return matrix; centered and standardised data.
#' @export
center_and_standardise <- function(my_data) {
  apply(my_data, MARGIN = 2, function(c1) {(c1 - mean(c1)) / sd(c1)})
}

# Used by simulate_data.R
rand <- function(nrow, ncol) {
  matrix(runif(nrow * ncol), ncol = ncol)
}
randn <- function(nrow, ncol) {
  matrix(rnorm(nrow * ncol), ncol = ncol)
}

L2_distance <- function(a, b, df) {
  if (nrow(a) == 1) {
    a = rbind(a, numeric(ncol(a)))
    b = rbind(b, numeric(ncol(b)))
  }
  aa = apply(a^2, 2, sum)  #vector
  bb = apply(b^2, 2, sum)  #vector
  ab = t(a) %*% b          #matrix

  m1 = matrix(rep(aa, length(bb)), ncol = length(bb))
  m2 = matrix(rep(bb, length(aa)), nrow = length(aa), byrow = TRUE)
  d = m1 + m2 - 2 * ab
  d = sqrt(d * (d > 0))
  if (df==1) {
    mmat = 1 - diag( max(nrow(d), ncol(d)) )
    d = d * ( mmat[1:nrow(d), 1:ncol(d)] );
  }
  d
}
