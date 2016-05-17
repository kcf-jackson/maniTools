library(pixmap)
library(magrittr)
load_files <- function(path) {
  file_list <- list.files(path)
  image_files <- purrr::map( file_list, ~ file.path(path, .) %>% read.pnm() )
  image_files
}
convert_files_to_matrix <- function(image_files) {
  purrr::map(image_files, getChannels)
}

## Utility functions
library(fields)
plot_faces <- function(i, proj_data, image_mat, imgw = 0.1, imgh = 0.1) {
  img <- image_mat[[i]]
  add.image(proj_data[i,1], proj_data[i,2], col = grey((0:256) / 256),
            t(img[nrow(img):1, ]), image.width = imgw, image.height = imgh)
}
