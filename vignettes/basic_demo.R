## ----setup, include=T----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = T, fig.width = 6, fig.height = 6)

## ------------------------------------------------------------------------
library(maniTools)
num_pts = 600
d = 2   #target dimension
k = 8   #k nearest neighbors
sim_data <- swiss_roll(num_pts)
plotly_3D(sim_data)

## ------------------------------------------------------------------------
# PCA (on centered and scaled data)
pca_dr <- sim_data$data %>% center_and_standardise() %>% prcomp()
proj_data <- sim_data$data %*% pca_dr$rotation[,1:2]
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "PCA")

# MDS
proj_data <- cmdscale(dist(sim_data$data), k = d)
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "MDS")

## ------------------------------------------------------------------------
# LLE
proj_data <- LLE2(sim_data$data, dim = d, k = k)
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "LLE")

# Hessian LLE
proj_data <- Hessian_LLE(sim_data$data, k = k, d = d)$projection
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "Hessian LLE")

# Laplacian Eigenmaps
proj_data <- Laplacian_Eigenmaps(sim_data$data, k = k, d = d)$eigenvectors
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "Laplacian Eigenmaps")

# LTSA
proj_data <- Local_TSA(sim_data$data, k = k, d = d)
plotly_2D(proj_data, sim_data$colors) %>% plotly::layout(title = "Local Tangent Space Alignment")

