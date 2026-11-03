#' Doubly Robust Win Statistics for Randomized Trials
#'
#' @name DRwinStats-package
#' @rdname DRwinStats-package
#' @aliases DRwinStats
#' @description Formula-based estimation of fixed-horizon prioritized win statistics
#' in individually randomized and cluster-randomized trials, with longitudinal
#' event records, flexible working survival models, and analytic variance estimation.
#' @details Start with \code{\link{drwins}} and \code{\link{drwins-data-format}}.
#' The introductory vignette is available through \code{vignette("DRwinStats")}.
#' Synthetic data are provided in \code{\link{example_irt}} and \code{\link{example_crt}}.
#' @author Xi Fang \email{xfang@mcw.edu}
#' @keywords package
#' @useDynLib DRwinStats, .registration = TRUE
#' @importFrom stats coef confint vcov setNames
"_PACKAGE"
