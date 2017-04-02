#' Swiss Roll
#' @param N integer; Number of datapoints.
#' @param height numeric; height of the Swiss Roll.
#' @export
swiss_roll <- function(N, height = 1) {
  p = (3 * pi/2) * (1 + 2 * runif(N, 0, 1))
  y = 21 * runif(N, 0, 1)
  ret_X = cbind(p * cos(p), y, height * p * sin(p))
  ret_ColorVector = p
  list(data = ret_X, colors = ret_ColorVector)
}



#' Swiss Hole
#' @param N integer; Number of datapoints.
#' @param height numeric; height of the Swiss Roll.
#' @export
swiss_hole <- function(N, height = 1) {
  p = (3 * pi/2) * (1 + 2 * runif(2*N, 0, 1))
  y = 21 * runif(2*N, 0, 1)
  kl = rep(0, 2 * N)
  assign_one = (p > 9) & (p < 12) & (y > 9) & (y < 14)
  kl[assign_one] = 1
  kkz = which(!assign_one)
  p = p[kkz[1:N]]
  y = y[kkz[1:N]]
  ret_X = cbind(p * cos(p), y, height * p * sin(p))
  ret_ColorVector = p
  list(data = ret_X, colors = ret_ColorVector)
}



#' Corner Planes
#' @param N integer; Number of datapoints.
#' @param angles numeric; Angle between the two planes.
#' @export
corner_planes <- function(N, angles = 45) {
  k = 1
  xMax = floor(sqrt(N))
  yMax = ceiling(N / xMax)
  cornerPoint = floor(yMax / 2)
  num_row = (xMax + 1) * (yMax + 1)
  X = matrix(numeric(num_row * 3), ncol = 3)
  ColorVector = numeric(num_row)
  for (x in 0:xMax) {
    for (y in 0:yMax) {
      if (y <= cornerPoint) {
        X[k, ] = c(x, y, 0)
      } else {
        newy = cornerPoint + (y - cornerPoint) * cos(pi * angles / 180)
        newz = (y - cornerPoint) * sin(pi * angles / 180)
        X[k, ] = c(x, newy, newz)
      }
      ColorVector[k] = y
      k = k+1;
    }
  }
  list(data = X, colors = ColorVector)
}



#' Punctured Sphere by Saul & Roweis
#' @param N integer; Number of datapoints.
#' @param z_scale numeric; "height" of the sphere.
#' @export
punctured_sphere <- function(N, z_scale = 1.0) {
  inc = 9 / sqrt(N);
  s = seq(-5, 5, by = inc)
  yy = matrix(rep(s, length(s)), nrow = length(s))
  xx = t(yy)
  rr2 = as.vector(xx)^2 + as.vector(yy)^2
  ii = order(rr2)
  Y = rbind(xx[ii[1:N]], yy[ii[1:N]])
  a = 4 / ( 4 + apply(Y^2, MARGIN = 2, sum) )
  X = cbind(a * Y[1,], a * Y[2,], z_scale * 2 * (1 - a))
  list(data = X, colors = X[,3])
}



#' Twin Peaks by Saul & Roweis
#' @param N integer; Number of datapoints.
#' @param z_scale numeric; "height" of the peaks.
#' @export
twin_peaks <- function(N, z_scale = 1.0) {
  inc = 1.5 / sqrt(N);
  s = seq(-1, 1, by = inc)
  yy2 = matrix(rep(s, length(s)), nrow = length(s))
  xx2 = t(yy2)
  zz2 = sin(pi * xx2) * tanh(3 * yy2)
  xy = 1 - 2 * rand(N, 2);
  X = cbind(xy, sin(pi * xy[,1] ) * tanh(3 * xy[,2]) * z_scale )
  list(data = X, colors = X[,3])
}



#' 3D Clusters
#' @param N integer; Number of datapoints.
#' @param num_cluster integer; number of clusters.
#' @export
clusters_3d <- function(N, num_cluster = 3) {
  numClusters = max(1, num_cluster)
  Centers = 10 * rand(numClusters, 3)

  D = L2_distance(t(Centers), t(Centers), 1)
  minDistance = min( D[D>0] )
  k = 1
  N2 = N - (numClusters - 1) * 9
  ColorVector = c()
  num_row = numClusters * ceiling(N2 / numClusters) + 9 * (num_cluster - 1)
  X = matrix(numeric(num_row * 3), ncol = 3)
  for (i in 1:numClusters) {
    for (j in 1:ceiling(N2 / numClusters)) {
      X[k, 1:3] = Centers[i,1:3] + (rand(1,3) - 0.5) * minDistance / sqrt(12)
      ColorVector = c(ColorVector, i)
      k = k + 1
    }
    # Connect clusters with straight line
    if (i < numClusters){
      for (t in seq(0.1, 0.9, by = 0.1)) {
        X[k,1:3] = Centers[i,1:3] + (Centers[i+1, 1:3] - Centers[i,1:3]) * t
        ColorVector = c(ColorVector, 0)
        k = k + 1;
      }
    }
  }
  list(data = X, colors = ColorVector)
}



#' Toroidal Helix by Coifman & Lafon
#' @param N integer; Number of datapoints.
#' @param sample_rate numeric; Sampling rate.
#' @export
toroidal_helix <- function(N, sample_rate = 1.0) {
  noiseSigma = 0.05 #noise parameter
  t = (1:N) / N
  t = t^(sample_rate) * 2 * pi
  X = cbind( (2+cos(8*t))*cos(t), (2+cos(8*t))*sin(t), sin(8*t) ) +
    noiseSigma * randn(N,3);
  list(data = X, colors = t)
}



#' Randomly sample from Gaussian distribution
#' @param N integer; Number of datapoints.
#' @param sigma numeric; Controls the variance of the normal distribution.
#' @export
gaussian_random_samples <- function(N, sigma = 1.0) {
  X = randn(N, 3)
  X[,3] = 1 / (sigma^2 * 2 * pi) * exp( (-X[,1]^2 - X[,2]^2) / (2 * sigma^2) )
  list(data = X, colors = X[,3])
}



#' A simple 3d plot function using plotly.
#' @param sim_data The data to be plotted.
#' @examples
#' sim_data <- swiss_roll(600)
#' library(plotly)
#' plotly_3D(sim_data)
#' @export
plotly_3D <- function(sim_data) {
  x <- sim_data$data[,1]
  y <- sim_data$data[,2]
  z <- sim_data$data[,3]
  scale <- sim_data$colors
  plotly::plot_ly(x = x, y = y, z = z, type = "scatter3d", mode = "markers",
          color = scale, marker = list(size = 5))
}



#' A simple 2d plot function using plotly.
#' @param proj_data The data to be plotted.
#' @param colors Color spectrum in real values. See example for usage.
#' @examples
#' sim_data <- swiss_roll(600)
#' library(plotly)
#' lle_data <- LLE2(sim_data$data, dim = 2, k = 8)
#' plotly_2D(lle_data, colors = sim_data$colors) #colors are optional
#' @export
plotly_2D <- function(proj_data, colors) {
  x <- proj_data[,1]
  y <- proj_data[,2]
  if (missing(colors))
    return(plotly::plot_ly(x = x, y = y, mode = "markers", marker = list(size = 5)))
  plotly::plot_ly(x = x, y = y, mode = "markers", color = colors, marker = list(size = 5))
}
