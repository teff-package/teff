#' teff (treatment effects classification): Package for classifying individuals into
#' subpoulations associated with positive and negative treatment effects using feature data.
#'
#' \code{teff} is a package designed to extract the profiles of subpopulations
#' with associated positive and negative treatment effect on an outcome variable. With the extracted profiles,
#' new individuals with feature data can be targeted and classified. If
#' treatment and outcome data is available for these new individuals the package can test
#' whether the association between the treatment and the outcome is indeed different
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
#' @import ggplot2
#' @import drc
#' @importFrom graphics lines points plot boxplot
#' @importFrom parallel mclapply
#' @importFrom betareg betareg
#' @importFrom stats lm glm as.formula complete.cases formula model.matrix predict
#' @importFrom survival coxph
#' @importFrom SummarizedExperiment assay
#' @importFrom ggpubr mean_ci ggline
#' @importFrom raster plot
#' @importFrom sva svaseq
#' @importFrom grDevices dev.off pdf
#' @importFrom limma voom
NULL
