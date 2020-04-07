#' Manifold learning in R
#' @description  This package is about porting Todd Wittman's Matlab codes, "MANI: Manifold Learning Toolkit", to R. In particular, three algorithms are ported: Laplacian Eigenmaps, Hessian LLE (HLLE) and Local Tangent Space Alignment (LTSA). Others are collected from existing R packages, including Isomap and LLE in "RDRToolbox", diffusionMap in "diffusionMap", t-SNE in "tsne", KPCA in "kernlab" and SPE in "spe".
#' @docType package
#' @name maniTools
#' @author Jackson Kwok
#' @references http://macs.citadel.edu/wittman/research.html
NULL
#> NULL

#' Pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Compound assignment pipe operator
#' @name %<>%
#' @rdname compound-assignment-pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

#' Import functions
#' @name imports
#' @keywords internal
#' @importFrom methods as
#' @importFrom stats rnorm runif sd
#' @importFrom utils tail
NULL
