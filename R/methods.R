.drw_conf_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("The confidence level must be strictly between zero and one.",
         call. = FALSE)
  as.numeric(level)
}

.drw_fitted_level <- function(object) {
  z <- object$settings$conf_level
  if (is.null(z)) z <- object$settings$conf.level
  if (is.null(z)) z <- .95
  .drw_conf_level(z)
}

.drw_measure_names <- function(measure) {
  if (is.null(measure)) return(NULL)
  if (!is.character(measure) || anyNA(measure))
    stop("measure must name a probability or win statistic.", call. = FALSE)
  key <- tolower(gsub("[ _-]", "", measure))
  lookup <- c(win = "win", loss = "loss", tie = "tie", nb = "NB",
              netbenefit = "NB", wr = "WR", winratio = "WR",
              wo = "WO", winodds = "WO")
  out <- unname(lookup[key])
  if (anyNA(out))
    stop("Unknown measure. Use win, loss, tie, NB, WR, or WO.",
         call. = FALSE)
  unique(out)
}

.drw_target_names <- function(target) {
  if (is.null(target)) return(NULL)
  if (!is.character(target) || anyNA(target))
    stop("target must be 'cluster' or 'individual'.", call. = FALSE)
  key <- tolower(gsub("[ _-]", "", target))
  lookup <- c(c = "cluster", cluster = "cluster", clusteraverage = "cluster",
              i = "individual", individual = "individual",
              individualaverage = "individual")
  out <- unname(lookup[key])
  if (anyNA(out))
    stop("target must be 'cluster' or 'individual'.", call. = FALSE)
  unique(out)
}

.drw_select_estimates <- function(object, measure = NULL,
                                   target = NULL, times = NULL) {
  out <- object$estimates
  if (!is.data.frame(out) ||
      !all(c("target", "time", "measure", "estimate") %in% names(out)))
    stop("The fitted object does not contain an estimates table.",
         call. = FALSE)
  measures <- .drw_measure_names(measure)
  targets <- .drw_target_names(target)
  if (!is.null(measures)) {
    if (length(setdiff(measures, as.character(out$measure))))
      stop("The requested measure was not estimated.", call. = FALSE)
    out <- out[as.character(out$measure) %in% measures, , drop = FALSE]
  }
  if (!is.null(targets)) {
    if (length(setdiff(targets, as.character(object$estimates$target))))
      stop("The requested target was not fitted. Available target(s): ",
           paste(unique(object$estimates$target), collapse = ", "), ".",
           call. = FALSE)
    out <- out[as.character(out$target) %in% targets, , drop = FALSE]
  }
  if (!is.null(times)) {
    if (!is.numeric(times) || anyNA(times) || any(!is.finite(times)) ||
        length(setdiff(times, object$estimates$time)))
      stop("times must be selected from the fitted restriction times.",
           call. = FALSE)
    out <- out[out$time %in% times, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.drw_interval_table <- function(object, out, level) {
  level <- .drw_conf_level(level)
  if (identical(level, .drw_fitted_level(object)) &&
      all(c("lower", "upper") %in% names(out))) return(out)
  z <- stats::qnorm((1 + level) / 2)
  if (!"se" %in% names(out)) out$se <- NA_real_
  out$lower <- out$upper <- rep(NA_real_, nrow(out))
  ratio <- as.character(out$measure) %in% c("WR", "WO")
  ok <- !ratio & is.finite(out$estimate) & is.finite(out$se) & out$se >= 0
  out$lower[ok] <- out$estimate[ok] - z * out$se[ok]
  out$upper[ok] <- out$estimate[ok] + z * out$se[ok]
  sel <- if ("se_log" %in% names(out)) out$se_log else
    ifelse(out$estimate > 0, out$se / out$estimate, NA_real_)
  ok <- ratio & is.finite(out$estimate) & out$estimate > 0 &
    is.finite(sel) & sel >= 0
  out$lower[ok] <- exp(log(out$estimate[ok]) - z * sel[ok])
  out$upper[ok] <- exp(log(out$estimate[ok]) + z * sel[ok])
  out
}

.drw_estimate_labels <- function(out) {
  paste(out$target, format(out$time, digits = 12, trim = TRUE),
        out$measure, sep = ":")
}

.drw_time_labels <- function(times, digits = 3L) {
  start <- max(1L, min(22L, as.integer(digits)))
  distinct <- unique(times)
  for (precision in seq.int(start, max(start, 12L))) {
    labels <- format(distinct, digits = precision, trim = TRUE)
    if (!anyDuplicated(labels))
      return(format(times, digits = precision, trim = TRUE))
  }
  sprintf("%.17g", times)
}

.drw_setting <- function(object, name, default = NULL) {
  out <- object$settings[[name]]
  if (is.null(out) && !is.null(object$call)) out <- object$call[[name]]
  if (is.null(out)) default else out
}

.drw_as_text <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.language(x) || inherits(x, "formula"))
    return(paste(deparse(x, width.cutoff = 80), collapse = " "))
  if (is.list(x)) {
    parts <- vapply(x, .drw_as_text, character(1))
    if (length(unique(parts)) == 1L) return(unname(parts[1L]))
    if (!is.null(names(x)) && all(nzchar(names(x))))
      parts <- paste(names(x), parts, sep = ": ")
    return(paste(parts, collapse = "; "))
  }
  paste(as.character(x), collapse = ", ")
}

.drw_model_convergence <- function(models) {
  codes <- numeric()
  flags <- logical()
  boundary <- FALSE
  visit <- function(z) {
    if (!is.list(z)) return(invisible(NULL))
    if (is.numeric(z$convergence))
      codes <<- c(codes, z$convergence)
    if (is.logical(z$converged))
      flags <<- c(flags, z$converged)
    if (isTRUE(z$boundary)) boundary <<- TRUE
    for (child in z[vapply(z, is.list, logical(1))]) visit(child)
    invisible(NULL)
  }
  visit(models)
  available <- length(codes) + length(flags) > 0L
  success <- available && !anyNA(codes) && !anyNA(flags) &&
    all(is.finite(codes)) && all(codes == 0) && all(flags)
  list(available = available, success = success, boundary = boundary)
}

.drw_print_header <- function(x, digits = 3L) {
  d <- x$design
  cat("Doubly robust win statistics\n")
  if (!is.null(d$type)) {
    if (d$type == "IRT") cat("Design: IRT; ", d$n, " participants\n", sep = "")
    else cat("Design: CRT; ", d$m, " clusters and ", d$n,
             " participants\n", sep = "")
  }
  if (!is.null(d$treatment_levels))
    cat("Comparison: ", d$treatment_levels[2L], " versus ",
        d$treatment_levels[1L], "\n", sep = "")
  if (!is.null(d$allocation))
    cat("Allocation: ", if (d$allocation == "complete")
      "fixed treatment totals" else "independent assignment", "\n", sep = "")
  if (!is.null(d$type) && d$type == "CRT" &&
      !is.null(d$cluster_size_summary)) {
    z <- d$cluster_size_summary
    cat("Cluster size: mean ", format(z["mean"], digits = digits, trim = TRUE),
        ", range ", z["min"], " to ", z["max"], "\n", sep = "")
  }
  cat("Restriction times: ",
      paste(.drw_time_labels(sort(unique(x$estimates$time)), digits),
            collapse = ", "), "\n", sep = "")
  if (isTRUE(d$default_times))
    cat("Times selected from participant follow-up quantiles.\n")
  om <- .drw_setting(x, "outcome_model")
  cm <- .drw_setting(x, "censoring_model")
  cp <- .drw_setting(x, "copula")
  if (!is.null(om)) {
    cat("Outcome model: ", .drw_as_text(om), sep = "")
    if (!is.null(cp)) cat("; copula: ", .drw_as_text(cp), sep = "")
    cat("\n")
  }
  if (!is.null(cm)) cat("Censoring model: ", .drw_as_text(cm), "\n", sep = "")
  of <- .drw_setting(x, "outcome_formula")
  if (is.null(of)) of <- .drw_setting(x, "outcome")
  cf <- .drw_setting(x, "censoring_formula")
  if (is.null(cf)) cf <- .drw_setting(x, "censoring")
  if (!is.null(of)) cat("Outcome formula: ", .drw_as_text(of), "\n", sep = "")
  if (!is.null(cf)) cat("Censoring formula: ", .drw_as_text(cf), "\n", sep = "")
  ctrl <- x$settings$control
  if (is.list(ctrl)) {
    penalties <- unlist(ctrl[c("censor_ridge", "outcome_ridge",
                              "association_ridge")], use.names = TRUE)
    if (length(penalties)) {
      if (all(penalties == 0)) cat("Ridge penalties: off\n")
      else cat("Ridge penalties: ",
               paste(names(penalties), format(penalties, digits = digits),
                     sep = "=", collapse = ", "), "\n", sep = "")
    }
  }
  conv <- x$diagnostics$converged
  if (is.logical(conv) && length(conv) && !anyNA(conv)) {
    cat("Model convergence: ", if (all(conv)) "achieved" else
      "one or more fits require review", "\n", sep = "")
  } else {
    conv <- .drw_model_convergence(x$diagnostics$models)
    if (conv$available)
      cat("Model convergence: ", if (conv$success)
        "all checks passed" else "one or more fits require review", "\n", sep = "")
    if (conv$boundary)
      cat("An outcome model reached a numerical parameter boundary; see diagnostics.\n")
  }
}

.drw_print_estimates <- function(out, digits) {
  wanted <- intersect(c("target", "time", "measure", "estimate", "se",
                         "se_log", "lower", "upper", "p_value"), names(out))
  tab <- out[, wanted, drop = FALSE]
  tab$time <- .drw_time_labels(tab$time, digits)
  if ("se_log" %in% names(tab) && all(is.na(tab$se_log)))
    tab$se_log <- NULL
  if ("p_value" %in% names(tab))
    tab$p_value[!tab$measure %in% c("NB", "WR", "WO")] <- NA_real_
  optional <- intersect(c("se", "se_log", "lower", "upper", "p_value"),
                        names(tab))
  for (nm in optional) if (all(is.na(tab[[nm]]))) tab[[nm]] <- NULL
  print(tab, row.names = FALSE, digits = digits, na.print = "")
}

#' Print and Extract Results from a DRwinStats Fit
#'
#' @name drwins-methods
#' @rdname drwins-methods
#' @aliases print.drwins coef.drwins vcov.drwins confint.drwins as.data.frame.drwins
#' @description Prints a concise description of a fitted \code{\link{drwins}} object or
#' extracts its estimates, confidence intervals, and covariance matrix.
#' The \code{\link{summary.drwins}} and \code{\link{plot.drwins}} help pages
#' describe the full summary and graphical display methods separately.
#' @param x,object A fitted object of class \code{"drwins"}, returned by
#' \code{\link{drwins}}. Use the generic functions \code{print}, \code{coef},
#' \code{vcov}, \code{confint}, and \code{as.data.frame} to call these methods.
#' @param digits Number of significant digits used by the print method. Time
#' labels use additional digits when needed to distinguish the fitted times.
#' The stored numerical values are unchanged.
#' @param measure A character vector containing one or more of \code{"win"},
#' \code{"loss"}, \code{"tie"}, \code{"NB"}, \code{"WR"}, or \code{"WO"}.
#' The default \code{NULL} selects every fitted measure. Names are
#' case-insensitive, and the full names of the three win statistics are also
#' accepted.
#' @param target \code{"cluster"}, \code{"individual"}, a vector containing
#' both, or \code{NULL} for all fitted targets. The selection must refer to
#' targets stored in the fitted object.
#' @param times A numeric vector selecting times from \code{object$times} or
#' \code{x$times}. The default \code{NULL} selects all fitted restriction times.
#' The methods return the existing results in their stored order.
#' @param parm For \code{confint}, an optional character vector of estimate
#' names or a numeric vector of row positions after selection by \code{measure},
#' \code{target}, and \code{times}. The estimate names have the same form as
#' the names returned by \code{coef}. The default \code{NULL} includes all
#' selected estimates.
#' @param level A single confidence level strictly between zero and one.
#' The default for \code{confint} is 0.95. The limits are calculated using
#' the stored estimates and standard errors.
#' @param row.names Optional row names for the returned data frame.
#' @param optional Included for compatibility with the \code{as.data.frame}
#' generic and currently unused.
#' @param ... Additional arguments for generic-method compatibility. They are
#' currently unused by these methods.
#' @details The print method displays the design, treatment comparison, model choices,
#' win ratio, and net benefit. It also reports overall censoring before the
#' largest fitted restriction time. Use \code{summary(fit)} for all six measures
#' and censoring rates by time, treatment arm, and target population.
#'
#' The \code{coef} method returns estimated probabilities or win statistics,
#' rather than the coefficients of the working outcome and censoring models.
#' Its names identify target, time, and measure, separated by colons, for example
#' \code{"individual:1:WR"}.
#'
#' The matrix returned by \code{vcov} is the covariance matrix of the estimated
#' win and loss probabilities, jointly across every fitted target and time.
#' It includes the sample-size scaling and the contribution from estimation of
#' the working models. Its dimensions do not correspond to all six measures
#' returned by \code{coef}. Row and column names identify the probability
#' coordinates. For standard errors of a selected win statistic, use the
#' \code{se} and \code{se_log} columns of \code{summary(fit)$estimates}.
#'
#' The \code{confint} method uses normal Wald intervals on the original scale
#' for probabilities and net benefit. It constructs win ratio and win odds
#' intervals on the logarithmic scale and exponentiates the limits. The
#' intervals are pointwise and are not adjusted for multiple measures or times.
#' The probability intervals are not clipped to zero or one.
#'
#' The \code{as.data.frame} method returns selected rows of the stored estimates
#' table, retaining the original confidence level. To obtain a table at a
#' different level, use \code{summary(fit, conf.level = 0.90)$estimates}.
#'
#' Fits with \code{inference = "none"} support printing, coefficient extraction,
#' and data-frame extraction. The \code{vcov} and \code{confint} methods require
#' an object fitted with \code{inference = "sandwich"}.
#' @return \describe{
#' \item{print}{Returns the fitted object invisibly after displaying it.}
#' \item{coef}{Returns a named numeric vector of selected estimated
#' probabilities and win statistics.}
#' \item{vcov}{Returns a numeric covariance matrix with one row and column for
#' each estimated win or loss probability at each fitted target and time.}
#' \item{confint}{Returns a two-column matrix of lower and upper confidence
#' limits, with estimate labels as row names.}
#' \item{as.data.frame}{Returns a data frame with columns \code{target},
#' \code{time}, \code{measure}, \code{estimate}, \code{se}, \code{se_log},
#' \code{lower}, \code{upper}, \code{p_value}, and \code{scale}. The
#' \code{estimate} and \code{se} columns use the original scale. The
#' \code{scale} column identifies the scale used to construct confidence
#' intervals.}
#' }
#' @seealso \code{\link{drwins}}, \code{\link{summary.drwins}},
#' \code{\link{plot.drwins}}, \code{\link{drwins_control}}
#' @examples data(example_irt, package = "DRwinStats")
#' final <- example_irt[example_irt$status %in% c(0L, 1L), ]
#' ids <- unlist(lapply(split(final$id, final$treatment),
#'                      function(z) head(z, 20L)), use.names = FALSE)
#' small <- example_irt[example_irt$id %in% ids, ]
#' fit <- drwins(
#'   small, outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km", copula = "independence",
#'   priorities = c(1L, 2L), times = c(1, 2), allocation = "complete"
#' )
#'
#' print(fit, digits = 4)
#' coef(fit, measure = c("WR", "NB"))
#' confint(fit, measure = "WR", level = 0.90)
#' vcov(fit)
#' as.data.frame(fit, measure = c("win", "loss", "tie"), times = 2)
#'
#' ## Select a confidence interval by its coefficient label.
#' label <- names(coef(fit, measure = "WR", times = 2))
#' confint(fit, parm = label)
#' @keywords models
#' @export
print.drwins <- function(x, digits = 3, ...) {
  .drw_print_header(x, digits)
  cat("\n")
  keep <- intersect(c("WR", "NB"), unique(as.character(x$estimates$measure)))
  .drw_print_estimates(
    .drw_select_estimates(x, if (length(keep)) keep else NULL), digits)
  cens <- x$censoring
  if (is.data.frame(cens) && nrow(cens)) {
    top <- max(cens$time)
    use <- cens$time == top & cens$arm == "Overall"
    if (any(use)) {
      top_label <- .drw_time_labels(cens$time, digits)[match(top, cens$time)]
      cat("\nCensored before time ", top_label, ": ",
          paste0(cens$target[use], " ",
                 formatC(100 * cens$rate[use], digits = 1, format = "f"), "%",
                 collapse = "; "), ".\n", sep = "")
    }
  }
  cat("Use summary() for all estimated probabilities and win statistics.\n")
  invisible(x)
}

#' Summarize Estimated Win Probabilities and Win Statistics
#'
#' @name summary.drwins
#' @rdname summary.drwins
#' @aliases print.summary.drwins
#' @description Summarizes a fitted \code{\link{drwins}} object for selected measures, target
#' populations, and restriction times. The summary reports estimated win,
#' loss, and tie probabilities, net benefit, win ratio, and win odds, together
#' with standard errors, confidence intervals, and censoring rates. Estimates
#' and censoring summaries are also returned as data frames for further use.
#' @param object A fitted object of class \code{"drwins"}, returned by
#' \code{\link{drwins}}. Call the method through \code{summary(object)}.
#' @param measure A character vector specifying the quantities to report.
#' The choices are \code{"win"}, \code{"loss"}, \code{"tie"}, \code{"NB"}
#' (net benefit), \code{"WR"} (win ratio), and \code{"WO"} (win odds).
#' The default \code{NULL} selects all six measures. Names are case-insensitive;
#' \code{"net benefit"}, \code{"win ratio"}, and \code{"win odds"} are also
#' accepted. A vector such as \code{c("WR", "NB")} selects several measures.
#' @param target A character vector containing \code{"cluster"},
#' \code{"individual"}, or both. The default \code{NULL} selects every target
#' stored in \code{object}. The cluster-average target gives equal weight to
#' clusters, and the individual-average target gives equal weight to participants.
#' An IRT fitted with the default \code{target = "both"} in \code{drwins} stores
#' the single individual-average target. Only fitted targets can be selected.
#' @param times A numeric vector selecting restriction times already present in
#' \code{object$times}. The default \code{NULL} retains all fitted times. Use the
#' stored numeric values, for example \code{object$times[1L]}. The method selects
#' existing results without interpolation or fitting the estimator at new times.
#' The selected rows retain their order in the fitted object.
#' @param conf.level A single confidence level strictly between zero and one.
#' The default \code{NULL} uses the level stored in the fit, which is 0.95 unless
#' changed through \code{conf_level} in \code{drwins}. For example,
#' \code{conf.level = 0.90} recalculates 90\% confidence intervals using the
#' stored estimates and standard errors. The argument to this summary method
#' is named \code{conf.level}, with a period.
#' @param digits Number of significant digits used when printing estimates and
#' censoring summaries. The default is 3. Additional digits are used for time
#' labels when needed to distinguish the fitted times. The stored numerical
#' results are not rounded.
#' @param x An object of class \code{"summary.drwins"}, returned by this summary
#' method.
#' @param ... Additional arguments for compatibility with the \code{summary}
#' and \code{print} generics. They are currently unused by these methods.
#' @details \subsection{Measures and interpretation}{
#' The treatment comparison follows the order specified by
#' \code{treatment_levels} in the fitted model: the second treatment level is
#' compared with the first. Write \eqn{p_W}, \eqn{p_L}, and \eqn{p_T} for the win,
#' loss, and tie probabilities at a particular restriction time and for a
#' particular target population. The reported summaries are
#' \eqn{\mathrm{NB} = p_W-p_L}, \eqn{\mathrm{WR} = p_W/p_L}, and
#' \eqn{\mathrm{WO} = (p_W+p_T/2)/(p_L+p_T/2)}.
#' The tie probability is calculated as \eqn{p_T = 1-p_W-p_L}.
#'
#' A positive net benefit or a win ratio or win odds above one favors the second
#' treatment level. The net benefit is an absolute difference between win and
#' loss probabilities. A win ratio compares those probabilities on a relative
#' scale. These quantities describe the prespecified hierarchical comparison.
#' Their definitions, clinical thresholds, and censoring adjustment are given
#' in \code{\link{drwins}}.
#' }
#' \subsection{Standard errors and confidence intervals}{
#' The \code{estimate} and \code{se} columns use the original scale of each
#' quantity. For the win ratio and win odds, \code{se_log} contains the standard
#' error of the logarithm. Their confidence limits are calculated on the
#' logarithmic scale and exponentiated. The \code{scale} column identifies the
#' scale used to construct the interval; it does not transform the values in
#' \code{estimate} or \code{se}.
#'
#' Probabilities and net benefit use normal Wald intervals on their original
#' scales. The \code{p_value} column gives two-sided normal tests of zero net
#' benefit, unit win ratio, or unit win odds, as applicable. No test is reported
#' for an individual win, loss, or tie probability. The intervals are pointwise
#' for each reported quantity and restriction time. Changing \code{conf.level}
#' uses the existing standard errors and leaves the fitted models and point
#' estimates unchanged.
#'
#' The summary retains the raw augmented estimates. Probability estimates and
#' their intervals are not constrained to the interval from zero to one. A
#' ratio or logarithmic interval that is undefined is represented by a missing
#' value. If the model was fitted with \code{inference = "none"}, the printed
#' summary reports point estimates and censoring rates. Inferential columns in
#' the returned estimates table contain missing values.
#' }
#' \subsection{Censoring rates}{
#' Censoring is reported before each selected restriction time, for each
#' treatment arm and overall. A participant is counted as censored at time
#' \code{tau} when the final record has status 0 and its time is strictly less
#' than \code{tau}. A status-0 record at or after \code{tau}, or an observed
#' terminal event, resolves the restricted outcome at that time.
#'
#' For the individual-average target, the rate is the proportion of participants
#' censored. For the cluster-average target, it is the mean of the within-cluster
#' censoring proportions, giving equal weight to clusters. The displayed
#' participant counts therefore need not give the cluster-average rate when
#' divided. The returned \code{rate} is a fraction, while the printed summary
#' displays it as \code{censoring_percent}. Selection by \code{target} and
#' \code{times} also selects the corresponding censoring rows.
#' }
#' @return \code{summary} returns a list of class \code{"summary.drwins"} with the
#' following components:
#' \describe{
#' \item{call}{The original call to \code{drwins}.}
#' \item{estimates}{A data frame with one row per selected target, restriction
#' time, and measure. Its columns are \code{target}, \code{time}, \code{measure},
#' \code{estimate}, \code{se}, \code{se_log}, \code{lower}, \code{upper},
#' \code{p_value}, and \code{scale}. Confidence limits correspond to
#' \code{conf.level}.}
#' \item{censoring}{A data frame for the selected targets and times, with columns
#' \code{target}, \code{time}, \code{arm}, \code{n}, \code{clusters},
#' \code{censored}, and \code{rate}. The counts refer to participants and
#' clusters in the indicated treatment arm or overall.}
#' \item{design}{The original trial design information, including the numbers of
#' participants and independent randomization units, treatment labels, and
#' cluster sizes.}
#' \item{models}{The fitted outcome and censoring models, or \code{NULL} when
#' model storage was disabled through \code{drwins_control(keep_models = FALSE)}.}
#' \item{diagnostics}{The original model, censoring-weight, and outcome-support
#' diagnostics. These retain the information for the complete fitted object.}
#' \item{settings}{The model choices, formulas, comparison rule, and other
#' settings retained in the fitted object.}
#' \item{conf.level}{The confidence level used for this summary.}
#' \item{digits}{The requested default printing precision.}
#' }
#' The print method displays the trial and model information, selected results,
#' and censoring table, and returns its input invisibly. Access the numerical
#' tables with \code{s$estimates} and \code{s$censoring} after
#' \code{s <- summary(fit)}. Use \code{vcov(fit)} to obtain the joint covariance
#' matrix of the estimated win and loss probabilities.
#' @seealso \code{\link{drwins}}, \code{\link{plot.drwins}},
#' \code{\link{coef.drwins}}, \code{\link{confint.drwins}},
#' \code{\link{vcov.drwins}}, \code{\link{example_irt}},
#' \code{\link{example_crt}}
#' @examples ## A small IRT example with 20 participants from each arm.
#' data(example_irt, package = "DRwinStats")
#' final <- example_irt[example_irt$status %in% c(0L, 1L), ]
#' ids <- unlist(lapply(split(final$id, final$treatment),
#'                      function(z) head(z, 20L)), use.names = FALSE)
#' small_irt <- example_irt[example_irt$id %in% ids, ]
#' fit_irt <- drwins(
#'   small_irt, outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence", priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 1.5, 2), allocation = "complete"
#' )
#'
#' ## All probabilities and win statistics, with censoring rates.
#' summary(fit_irt)
#'
#' ## Select measures and fitted times, and use 90% confidence intervals.
#' s <- summary(fit_irt, measure = c("WR", "NB"), times = c(1, 2),
#'              conf.level = 0.90)
#' print(s, digits = 4)
#' s$estimates
#' s$censoring
#'
#' ## A small CRT example retains four complete clusters from each arm.
#' data(example_crt, package = "DRwinStats")
#' cluster_arm <- unique(example_crt[c("cluster", "treatment")])
#' keep <- unlist(lapply(split(cluster_arm$cluster, cluster_arm$treatment),
#'                       function(z) head(z, 4L)), use.names = FALSE)
#' small_crt <- example_crt[example_crt$cluster %in% keep, ]
#' fit_crt <- drwins(
#'   small_crt, cluster = "cluster", outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence", priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 2), target = "both", allocation = "complete"
#' )
#' summary(fit_crt, measure = c("WR", "WO", "NB"))
#' summary(fit_crt, measure = c("win", "loss", "tie"),
#'         target = "individual", times = 2)
#' @keywords models
#' @export
summary.drwins <- function(object, measure = NULL, target = NULL, times = NULL,
                           conf.level = NULL, digits = 3, ...) {
  if (is.null(conf.level)) conf.level <- .drw_fitted_level(object)
  estimates <- .drw_interval_table(
    object, .drw_select_estimates(object, measure, target, times), conf.level)
  censoring <- object$censoring
  if (is.data.frame(censoring)) {
    censoring <- censoring[censoring$target %in% unique(estimates$target) &
                             censoring$time %in% unique(estimates$time),
                           , drop = FALSE]
    rownames(censoring) <- NULL
  }
  ans <- list(call = object$call, estimates = estimates, censoring = censoring,
              design = object$design, models = object$models,
              diagnostics = object$diagnostics, settings = object$settings,
              conf.level = .drw_conf_level(conf.level), digits = digits)
  class(ans) <- "summary.drwins"
  ans
}

#' @rdname summary.drwins
#' @export
print.summary.drwins <- function(x, digits = x$digits, ...) {
  .drw_print_header(x, digits)
  if (identical(x$settings$inference, "none"))
    cat("\nPoint estimates (inference was not requested)\n")
  else
    cat("\nEstimates and ", format(100 * x$conf.level, trim = TRUE),
        "% confidence intervals\n", sep = "")
  .drw_print_estimates(x$estimates, digits)
  if (is.data.frame(x$censoring) && nrow(x$censoring)) {
    cat("\nCensoring before each restriction time\n")
    z <- x$censoring
    z$time <- .drw_time_labels(z$time, digits)
    z$censoring_percent <- 100 * z$rate
    z$rate <- NULL
    print(z, row.names = FALSE, digits = digits)
    if ("cluster" %in% z$target)
      cat("Cluster-target rates average the participant proportions within clusters.\n")
    cat("A follow-up ending exactly at the restriction time is resolved at that time.\n")
  }
  invisible(x)
}

#' @rdname drwins-methods
#' @export
coef.drwins <- function(object, measure = NULL, target = NULL, times = NULL, ...) {
  z <- .drw_select_estimates(object, measure, target, times)
  stats::setNames(z$estimate, .drw_estimate_labels(z))
}

#' @rdname drwins-methods
#' @export
vcov.drwins <- function(object, ...) {
  if (identical(object$settings$inference, "none") ||
      is.null(object$covariance))
    stop("The fit does not contain a covariance matrix; fit with inference = 'sandwich'.",
         call. = FALSE)
  object$covariance
}

#' @rdname drwins-methods
#' @export
confint.drwins <- function(object, parm = NULL, level = .95, measure = NULL,
                           target = NULL, times = NULL, ...) {
  if (identical(object$settings$inference, "none"))
    stop("Confidence intervals were not calculated; fit with inference = 'sandwich'.",
         call. = FALSE)
  z <- .drw_interval_table(
    object, .drw_select_estimates(object, measure, target, times), level)
  labels <- .drw_estimate_labels(z)
  if (!is.null(parm)) {
    if (is.character(parm)) {
      idx <- match(parm, labels)
      if (anyNA(idx))
        stop("Unknown parm; use names(coef(object)) to inspect estimate names.",
             call. = FALSE)
    } else if (is.numeric(parm) && !anyNA(parm) &&
               all(parm == floor(parm) & parm >= 1 & parm <= nrow(z))) {
      idx <- as.integer(parm)
    } else stop("parm must select estimate names or row positions.", call. = FALSE)
    z <- z[idx, , drop = FALSE]
    labels <- labels[idx]
  }
  out <- as.matrix(z[, c("lower", "upper"), drop = FALSE])
  rownames(out) <- labels
  colnames(out) <- paste0(format(100 * c((1 - level) / 2,
                                       (1 + level) / 2), trim = TRUE), " %")
  out
}

#' @rdname drwins-methods
#' @export
as.data.frame.drwins <- function(x, row.names = NULL, optional = FALSE,
                                 measure = NULL, target = NULL, times = NULL, ...) {
  out <- .drw_select_estimates(x, measure, target, times)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' Plot Win Probabilities or Win Statistics at Fitted Restriction Times
#'
#' @name plot.drwins
#' @rdname plot.drwins
#' @description Creates a base R plot of a probability or win statistic from a fitted
#' \code{\link{drwins}} object. Estimates are displayed at the fitted restriction
#' times, with optional pointwise confidence intervals. A single display can
#' compare the cluster-average and individual-average targets, or show all
#' three comparison probabilities.
#' @param x A fitted object of class \code{"drwins"}, returned by
#' \code{\link{drwins}}. Call the method through \code{plot(x)}.
#' @param measure A single character string naming the quantity to display:
#' \code{"win"}, \code{"loss"}, \code{"tie"}, \code{"NB"} (net benefit),
#' \code{"WR"} (win ratio), or \code{"WO"} (win odds). The default is
#' \code{"WR"}. Use \code{"probabilities"} to display win, loss, and tie
#' probabilities together. The individual measure names are case-insensitive,
#' and the full names \code{"net benefit"}, \code{"win ratio"}, and
#' \code{"win odds"} are also accepted.
#' @param target A character vector containing \code{"cluster"},
#' \code{"individual"}, or both. The default \code{NULL} displays every target
#' stored in \code{x}. Only fitted targets can be selected. An IRT fitted using
#' the default settings has the single individual-average target.
#' @param conf.int Logical. If \code{TRUE}, the default, draw confidence
#' intervals where finite interval limits are available. Set to \code{FALSE}
#' to display the estimates alone. A fit with \code{inference = "none"}
#' provides point estimates without interval bars.
#' @param level A single confidence level strictly between zero and one.
#' The default \code{NULL} uses the level stored in the fit, which is 0.95 unless
#' \code{conf_level} was changed in \code{drwins}. For example,
#' \code{level = 0.90} displays 90\% confidence intervals using the existing
#' standard errors. The argument to \code{summary.drwins} is instead named
#' \code{conf.level}.
#' @param xlab Label for the horizontal axis. The default is
#' \code{"Restriction time"}. Times have the same units as the input event
#' and follow-up times.
#' @param ylab Label for the vertical axis. The default \code{NULL} supplies
#' \code{"Win ratio"}, \code{"Win odds"}, \code{"Net benefit"},
#' \code{"Win probability"}, \code{"Loss probability"}, or
#' \code{"Tie probability"}, according to the selected measure. With
#' \code{measure = "probabilities"}, it supplies \code{"Comparison probability"}.
#' @param col Colors for the plotted groups. A group is one combination of a
#' target population and a measure. Values are recycled in the order of groups
#' in the selected estimates table. By default, win, loss, and tie probabilities
#' use muted sage, terracotta, and gray, respectively. Other displays use a
#' muted palette to distinguish targets. The interval bars use a translucent
#' version of each group's color.
#' @param pch Point symbols, using the conventions of
#' \code{\link[graphics]{points}}. Values are recycled across groups. The
#' defaults distinguish target populations with different symbols.
#' @param lty Line types, using the conventions of
#' \code{\link[graphics]{lines}}. Values are recycled across groups. The
#' defaults distinguish target populations with different line types.
#' @param type The default \code{"b"} draws both points and connecting lines.
#' Use \code{"p"} for points alone or \code{"l"} for lines alone.
#' @param legend Logical. The default \code{TRUE} adds a legend in the upper
#' right when more than one group is displayed, including displays of all
#' three probabilities. Set to \code{FALSE} to omit it or add a separate legend
#' with \code{\link[graphics]{legend}}.
#' @param log Axes to display logarithmically, as in
#' \code{\link[graphics]{plot.default}}: \code{""} for ordinary axes,
#' \code{"x"} for the horizontal axis, \code{"y"} for the vertical axis,
#' or \code{"xy"} for both. The default is \code{""}, including for ratio
#' measures. A logarithmic vertical axis requires positive plotted estimates.
#' @param ... Additional arguments passed to the initial
#' \code{\link[graphics]{plot.default}} call. Examples include \code{xlim},
#' \code{ylim}, \code{main}, \code{cex.axis}, and \code{cex.lab}. Use the named
#' \code{col}, \code{pch}, and \code{lty} arguments to control the plotted
#' series. The additional arguments are not forwarded to the later calls that
#' draw individual points, lines, or interval bars.
#' @details \subsection{Comparison and target populations}{
#' The treatment comparison is the second treatment level versus the first,
#' as recorded in \code{x$design$treatment_levels}. Positive net benefit and
#' win ratio or win odds above one favor the second treatment level. A horizontal
#' reference line is added at one for the win ratio and win odds and at zero for
#' the net benefit when compatible with the selected axis scale.
#'
#' The cluster-average target gives equal weight to clusters, and the
#' individual-average target gives equal weight to participants. When both were
#' fitted, \code{target = NULL} displays both. With
#' \code{measure = "probabilities"}, color identifies the win, loss, or tie
#' probability by default, while point symbols and line types identify the
#' target population. See \code{\link{summary.drwins}} for definitions and
#' interpretation of the reported measures.
#' }
#' \subsection{Restriction times and intervals}{
#' All fitted restriction times are used. Connecting lines join the estimates
#' in increasing time order and do not add a smoothing model or estimate the
#' comparison at intervening times. To evaluate a different set of times, supply
#' those times to \code{drwins}. The graphical argument \code{xlim} changes the
#' visible range of the display without changing the fitted results.
#'
#' The displayed intervals are pointwise, using the analytic standard errors
#' stored in the fitted object. Win ratio and win odds intervals are calculated
#' on the logarithmic scale and exponentiated. Probability and net benefit
#' intervals use their original scales. Changing \code{level} recalculates the
#' limits from those standard errors without refitting the model. Setting
#' \code{log = "y"} changes the axis display; it does not change how the
#' estimator or confidence intervals are calculated.
#'
#' Probability estimates and their intervals are shown on their original
#' scale and are not clipped to zero or one. Undefined ratio results are omitted.
#' A nonpositive ratio estimate produces a warning and is omitted from the
#' display. On a logarithmic vertical axis, interval bars require positive
#' finite endpoints. The default limits include the finite displayed results
#' and, where applicable, the reference value. Supplying \code{ylim} controls
#' the visible range.
#' }
#' \subsection{Using the plotted values}{
#' The method uses base R graphics and returns the selected estimates table
#' invisibly. Assign the result, as in \code{z <- plot(fit, measure = "NB")},
#' to use the plotted values in another display or report. For a numerical
#' summary without opening a plot, use \code{\link{summary.drwins}} or
#' \code{\link{as.data.frame.drwins}}.
#' }
#' @return Invisibly returns a data frame containing the rows used for the display. Its
#' columns are \code{target}, \code{time}, \code{measure}, \code{estimate},
#' \code{se}, \code{se_log}, \code{lower}, \code{upper}, \code{p_value}, and
#' \code{scale}. The confidence limits correspond to the requested \code{level}.
#' Nonpositive ratio estimates omitted from the plot and interval endpoints
#' incompatible with a logarithmic axis are represented by missing values in
#' this returned table. The original fitted object is unchanged.
#' @seealso \code{\link{drwins}}, \code{\link{summary.drwins}},
#' \code{\link{as.data.frame.drwins}}, \code{\link[graphics]{plot.default}},
#' \code{\link{example_irt}}, \code{\link{example_crt}}
#' @examples ## A small IRT example with all observed records for 20 participants per arm.
#' data(example_irt, package = "DRwinStats")
#' final <- example_irt[example_irt$status %in% c(0L, 1L), ]
#' ids <- unlist(lapply(split(final$id, final$treatment),
#'                      function(z) head(z, 20L)), use.names = FALSE)
#' small_irt <- example_irt[example_irt$id %in% ids, ]
#' fit_irt <- drwins(
#'   small_irt, outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence", priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 1.5, 2), allocation = "complete"
#' )
#'
#' ## Default win ratio plot, including pointwise 95% confidence intervals.
#' plot(fit_irt)
#'
#' ## All comparison probabilities, with estimates alone.
#' plot(fit_irt, measure = "probabilities", conf.int = FALSE)
#'
#' ## Custom labels and styling, and a different confidence level.
#' plotted <- plot(fit_irt, measure = "NB", level = 0.90, type = "p",
#'                 col = "#A8654E", pch = 16,
#'                 xlab = "Years since randomization", ylab = "Net benefit")
#' plotted[c("time", "estimate", "lower", "upper")]
#'
#' ## Preserve whole clusters when constructing a small CRT example.
#' data(example_crt, package = "DRwinStats")
#' cluster_arm <- unique(example_crt[c("cluster", "treatment")])
#' keep <- unlist(lapply(split(cluster_arm$cluster, cluster_arm$treatment),
#'                       function(z) head(z, 4L)), use.names = FALSE)
#' small_crt <- example_crt[example_crt$cluster %in% keep, ]
#' fit_crt <- drwins(
#'   small_crt, cluster = "cluster", outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence", priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 2), target = "both", allocation = "complete"
#' )
#'
#' ## Compare the two fitted target populations.
#' plot(fit_crt, measure = "WR", col = c("#A8654E", "#5E8373"),
#'      pch = c(16, 17), lty = c(1, 2))
#'
#' ## A single target on a logarithmic ratio axis.
#' plot(fit_crt, measure = "WR", target = "individual", log = "y")
#' @keywords hplot
#' @export
plot.drwins <- function(x, measure = "WR", target = NULL, conf.int = TRUE,
                        level = NULL, xlab = "Restriction time", ylab = NULL,
                        col = NULL, pch = NULL, lty = NULL,
                        type = c("b", "p", "l"), legend = TRUE, log = "", ...) {
  type <- match.arg(type)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int) ||
      !is.logical(legend) || length(legend) != 1L || is.na(legend))
    stop("conf.int and legend must be TRUE or FALSE.", call. = FALSE)
  probabilities <- identical(tolower(measure), "probabilities")
  measures <- if (probabilities) c("win", "loss", "tie") else
    .drw_measure_names(measure)
  if (is.null(measures) || (!probabilities && length(measures) != 1L))
    stop("Choose one measure, or measure = 'probabilities'.", call. = FALSE)
  if (is.null(level)) level <- .drw_fitted_level(x)
  z <- .drw_interval_table(x, .drw_select_estimates(x, measures, target), level)
  if (!nrow(z)) stop("There are no selected estimates to plot.", call. = FALSE)
  if (!"lower" %in% names(z)) z$lower <- NA_real_
  if (!"upper" %in% names(z)) z$upper <- NA_real_
  ratio <- all(measures %in% c("WR", "WO"))
  if (ratio && any(is.finite(z$estimate) & z$estimate <= 0)) {
    warning("Nonpositive ratio estimates are not plotted; inspect summary(x).",
            call. = FALSE)
    bad <- is.finite(z$estimate) & z$estimate <= 0
    z$estimate[bad] <- z$lower[bad] <- z$upper[bad] <- NA_real_
  }
  if (grepl("y", log, fixed = TRUE)) {
    if (any(is.finite(z$estimate) & z$estimate <= 0))
      stop("A logarithmic y-axis requires positive plotted estimates.",
           call. = FALSE)
    z$lower[z$lower <= 0] <- NA_real_
    z$upper[z$upper <= 0] <- NA_real_
  }
  yrange <- z$estimate
  if (conf.int) yrange <- c(yrange, z$lower, z$upper)
  yrange <- yrange[is.finite(yrange)]
  if (!length(yrange)) stop("There are no finite estimates to plot.", call. = FALSE)
  reference <- if (ratio) 1 else if (identical(measures, "NB")) 0 else NA_real_
  if (is.finite(reference) &&
      !(grepl("y", log, fixed = TRUE) && reference <= 0))
    yrange <- c(yrange, reference)
  group_key <- paste(z$target, z$measure, sep = "\036")
  group <- factor(group_key, levels = unique(group_key))
  groups <- levels(group)
  ng <- length(groups)
  group_measures <- vapply(groups, function(g)
    as.character(z$measure[which(group == g)[1L]]), character(1))
  group_targets <- vapply(groups, function(g)
    as.character(z$target[which(group == g)[1L]]), character(1))
  palette <- c("#A8654E", "#5E8373", "#B29A59", "#6B7480", "#907364", "#4E716B")
  if (is.null(col)) col <- if (probabilities)
    unname(c(win = "#5E8373", loss = "#A8654E",
             tie = "#8A8B81")[group_measures]) else palette
  if (is.null(pch)) pch <- if (probabilities)
    c(16, 17)[match(group_targets, unique(group_targets))] else c(16, 17, 15, 18)
  if (is.null(lty)) lty <- if (probabilities)
    match(group_targets, unique(group_targets)) else c(1, 2, 3)
  col <- rep(col, length.out = ng)
  pch <- rep(pch, length.out = ng)
  lty <- rep(lty, length.out = ng)
  if (is.null(ylab)) ylab <- if (probabilities) "Comparison probability" else
    switch(measures, WR = "Win ratio", WO = "Win odds", NB = "Net benefit",
           win = "Win probability", loss = "Loss probability",
           tie = "Tie probability")
  dots <- list(...)
  if (is.null(dots$xlim)) {
    xr <- range(z$time)
    if (diff(xr) == 0) xr <- xr + c(-1, 1) * max(.05 * abs(xr[1]), .05)
    if (grepl("x", log, fixed = TRUE) && xr[1] <= 0)
      xr <- range(z$time) * c(.9, 1.1)
    dots$xlim <- xr
  }
  if (is.null(dots$ylim)) {
    yr <- range(yrange)
    if (diff(yr) == 0) {
      delta <- max(.05 * abs(yr[1]), .025)
      yr <- yr + c(-delta, delta)
      if (grepl("y", log, fixed = TRUE) && yr[1] <= 0)
        yr <- range(yrange) * c(.9, 1.1)
    }
    dots$ylim <- yr
  }
  do.call(graphics::plot.default,
          c(list(x = range(z$time), y = range(yrange), type = "n",
                 xlab = xlab, ylab = ylab, log = log), dots))
  if (is.finite(reference) &&
      !(grepl("y", log, fixed = TRUE) && reference <= 0))
    graphics::abline(h = reference, col = "#B7B5B0", lty = 3)
  legend_labels <- character(ng)
  for (j in seq_len(ng)) {
    use <- which(group == groups[j])
    use <- use[order(z$time[use])]
    zz <- z[use, , drop = FALSE]
    if (conf.int) {
      ci <- is.finite(zz$lower) & is.finite(zz$upper) & zz$upper > zz$lower
      if (any(ci))
        graphics::arrows(zz$time[ci], zz$lower[ci], zz$time[ci], zz$upper[ci],
                         length = .045, angle = 90, code = 3,
                         col = grDevices::adjustcolor(col[j], alpha.f = .65))
    }
    if (type %in% c("b", "l"))
      graphics::lines(zz$time, zz$estimate, col = col[j], lty = lty[j], lwd = 1.5)
    if (type %in% c("b", "p"))
      graphics::points(zz$time, zz$estimate, col = col[j], pch = pch[j])
    pop <- if (zz$target[1L] == "cluster") "Cluster average" else "Individual average"
    legend_labels[j] <- if (probabilities)
      paste(pop, zz$measure[1L], sep = ": ") else pop
  }
  if (legend && (ng > 1L || probabilities))
    graphics::legend("topright", legend = legend_labels, col = col,
                     lty = if (type == "p") NA else lty,
                     pch = if (type == "l") NA else pch,
                     lwd = 1.5, bty = "n", cex = .85)
  invisible(z)
}
