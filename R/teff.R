#' teff (treatment effects classification): Package for classifying individuals into
#' subpoulations associated with high and low treatment effects using feature data.
#'
#' \code{teff} is a package designed to extract the profiles of subpopulations
#' with associated high and low treatment effects. With the extracted profiles,
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
#' @importFrom parallel mclapply
#' @importFrom betareg betareg
#' @importFrom stats lm glm
#' @importFrom survival coxph
#' @importFrom SummarizedExperiment assay
#'
NULL
