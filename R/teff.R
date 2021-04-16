#' teff (treatment effects classification): Package for classifying individuals into
#' subpoulations associated with positive and negative treatment effects using feature data.
#'
#' \code{teff} is a package designed to extract the profiles of subpopulations
#' with associated positive and negative treatment effects. With the extracted profiles,
#' new individuals with feature data can be targeted and classified. If
#' treatment and effect data is available for these new individuals the package can test
#' whether the association between the treatment and the effect is indeed different
#' across subpopulations.
#'
#' @docType package
#' @name teff
#'
#' @import methods
#' @import grf
#' @import survival
#' @import Biobase
#' @import methods
#' @import graphics
#' @import minfi
#' @import ggplot2
#' @importFrom graphics lines points plot
#' @importFrom parallel mclapply
#' @importFrom betareg betareg
#' @importFrom stats lm glm
#' @importFrom survival coxph
#' @importFrom SummarizedExperiment assay
#' @importFrom ggpubr mean_ci ggline
#' @importFrom raster plot


NULL
