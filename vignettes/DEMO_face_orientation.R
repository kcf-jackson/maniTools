## ---- message = F, fig.width = 7, fig.height = 7-------------------------
rm(list = ls())
library(maniTools)
library(R.matlab)
util_path <- system.file("face_orientation", "util_face_orientation.R", package = "maniTools")
source(util_path)

# load data
fpath <- system.file("face_orientation", "face_data.mat", package = "maniTools")
matlab_file <- readMat(fpath)
image_mat <- t(matlab_file$images)
 
# perform dimensionality reduction
proj_data <- RDRToolbox::Isomap(image_mat, 2, 6)$dim2  #set d = 3 to learn lighting directions.
plot(proj_data, pch = 19, cex = 0.5)

# plot results
sample_image_ind <- sample(seq(nrow(image_mat)), 20)
purrr::walk(sample_image_ind, ~plot_faces(., proj_data, image_mat))

