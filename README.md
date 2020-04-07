## Manifold learning in R


### Installation

To install the package, run the following codes in R:
```
BiocManager::install("RDRToolbox")

# Without vignettes
devtools::install_github("kcf-jackson/maniTools")

# With vignettes (Note: it will take some time to build all the vignettes)
install.packages(c("knitr", "fields", "purrr", "animation", "R.matlab", "pixmap")
devtools::install_github("kcf-jackson/maniTools", build_vignettes = TRUE)
browseVignettes("maniTools")
```


### Pre-built vignettes

1. [A list of supported dimension-reduction techniques](https://raw.githack.com/kcf-jackson/maniTools/master/vignettes/basic_demo.html)
2. [Face ordering (Wang, 2011)](https://rawcdn.githack.com/kcf-jackson/maniTools/2eeb0ab18577dbd30c556390098d5c810da243dc/vignettes/DEMO_face_ordering.html))
3. [Face orientation (Tenenbaum, Silva and Langford, 2000)](https://rawcdn.githack.com/kcf-jackson/maniTools/2eeb0ab18577dbd30c556390098d5c810da243dc/vignettes/DEMO_face_orientation.html)
4. [Robot Wifi trajectory (Lawrence, 2012)](https://rawcdn.githack.com/kcf-jackson/maniTools/2eeb0ab18577dbd30c556390098d5c810da243dc/vignettes/DEMO_robot_wifi.html)
5. [Shiny App](https://raw.githack.com/kcf-jackson/maniTools/master/vignettes/shiny_app.html)


### References:

1. Wang, J. (2011). Geometric structure of high-dimensional data and dimensionality reduction (pp. 294). Springer Berlin Heidelberg.
2. Tenenbaum, J. B., De Silva, V., & Langford, J. C. (2000). A global geometric framework for nonlinear dimensionality reduction. science, 290(5500), 2319-2323.
3. Lawrence, N. D. (2012). A unifying probabilistic perspective for spectral dimensionality reduction: Insights and new models. The Journal of Machine Learning Research, 13(1), 1609-1638.
