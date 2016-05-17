## Utility functions
reshape_to_matrix <- function(vec, nr, nc, byrow = TRUE) {
  matrix(rev(vec), nrow = nr, ncol = nc, byrow = byrow)
}

library(fields)
plot_faces <- function(i, proj_data, image_mat) {
  add.image(proj_data[i,1], proj_data[i,2], col = grey((0:256) / 256),
            reshape_to_matrix(image_mat[i,], 64, 64, TRUE),
            image.width = 0.1, image.height = 0.1)
}
