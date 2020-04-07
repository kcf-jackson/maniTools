## ---- fig.height=3, fig.width=4, message=FALSE---------------------------
# Load data
rm(list = ls())
library(magrittr)
library(plotly)

fpath <- system.file("robot_wireless", "wireless_complete.csv", package = "maniTools")
wireless_full <- read.csv(fpath, header = FALSE)

ground_truth <- wireless_full[, 1:3]
names(ground_truth) <- c("x", "y", "time")

## ------------------------------------------------------------------------
fpath <- system.file("robot_wireless", "wireless.csv", package = "maniTools")
wireless <- read.csv(fpath, header = FALSE) %>% as.matrix()
dim(wireless)  # 215 datapoints of 30 wifi signal strength.

## ------------------------------------------------------------------------
# Perform dimensionality reduction
proj_data <- RDRToolbox::Isomap(wireless, dim = 2, k = 7)[[1]] %>% as.data.frame()
names(proj_data) <- c('x', 'y')

## ---- fig.height=4, fig.width=8------------------------------------------
# Plotting
p1 <- plot_ly(ground_truth[c('x', 'y')], x = ~x, y = ~y, name = "Truth")
p2 <- plot_ly(data = proj_data, x = ~x, y = ~y, mode = 'markers+lines', name = "Isomap")
subplot(p1, p2) %>% layout(title = "Isomap")

## ---- fig.height=4, fig.width=9, message = FALSE-------------------------
library(maniTools)
# Laplacian Eigenmaps
LE_proj_data <- Laplacian_Eigenmaps(wireless, k = 7, d = 2)[[1]] %>% as.data.frame()
names(LE_proj_data) <- c('x', 'y')

# Locally Linear Embedding
LLE_proj_data <- LLE2(wireless, dim = 2, k = 7) %>% as.data.frame()
names(LLE_proj_data) <- c('x', 'y')

p3 <- plot_ly(data = LE_proj_data, x = ~x, y = ~y, mode = 'markers+lines', name = "Laplacian Eigenmaps")
p4 <- plot_ly(data = LLE_proj_data, x = ~x, y = ~y, mode = 'markers+lines', name = "Locally Linear Embedding")
subplot(p3, p4)

