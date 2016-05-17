## ------------------------------------------------------------------------
# Load data
rm(list = ls())
library(magrittr)
fpath <- system.file("robot_wireless", "wireless_complete.csv", package = "maniTools")
wireless_full <- read.csv(fpath, header = FALSE)
ground_truth <- wireless_full[, 1:3]
names(ground_truth) <- c("x", "y", "time")

fpath <- system.file("robot_wireless", "wireless.csv", package = "maniTools")
wireless <- read.csv(fpath, header = FALSE) %>% as.matrix()
dim(wireless)  # 215 datapoints of 30 wifi signal strength.


# Perform dimensionality reduction
proj_data <- RDRToolbox::Isomap(wireless, dim = 2, k = 7)[[1]] %>% as.data.frame()
names(proj_data) <- c('x', 'y')
proj_data


# Plotting
library(plotly)
p1 <- plot_ly(ground_truth[c('x', 'y')], x = x, y = y, name = "Truth")
p2 <- plot_ly(data = proj_data, x = x, y = y, mode = 'markers + lines', name = "Model")
subplot(p1, p2) %>% layout(title = "Isomap")

