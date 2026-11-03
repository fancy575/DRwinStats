#' Control Fitting and Storage
#'
#' @name drwins_control
#' @rdname drwins_control
#' @description Options for model fitting, numerical accuracy, and storage
#' in \code{\link{drwins}}.
#' @param ridge Nonnegative default ridge penalty for outcome and censoring
#' regression coefficients on their original coefficient scale.
#' Zero requests unpenalized fitting.
#' @param censor_ridge Nonnegative ridge penalty for censoring regression
#' coefficients on their original coefficient scale.
#' @param outcome_ridge Nonnegative ridge penalty for outcome regression
#' coefficients on their original coefficient scale.
#' @param association_ridge Nonnegative penalty toward independence for
#' the association parameter in a dependent joint copula outcome fit.
#' @param baseline_penalty Nonnegative penalty on the final cumulative
#' baseline hazards in a dependent joint copula outcome fit.
#' @param maxit Positive integer maximum number of optimization iterations.
#' @param tol Positive convergence tolerance.
#' @param finite_sample Logical, with default \code{TRUE}, indicating whether
#' to apply finite sample adjustments to the empirical kernel projection and
#' covariance calculation. The adjustment uses the number of independent
#' participants or clusters, as described under \sQuote{Details}.
#' @param weight_floor Positive lower bound, less than one, used to check
#' fitted censoring survival probabilities.
#' @param positivity With \code{"error"}, a censoring survival probability
#' below \code{weight_floor} stops fitting. With \code{"truncate"}, the
#' probability is replaced by the floor and the fit records the truncation.
#' @param max_states Positive whole-number limit on the outcome states
#' used by the numerical calculation. It is at most
#' \code{.Machine$integer.max}. Exceeding the limit stops the calculation
#' rather than dropping an endpoint.
#' @param max_bytes Positive upper bound on the estimated memory requirement, in bytes.
#' @param keep_models Whether to retain fitted nuisance models in the result.
#' @param keep_influence Whether to retain the estimated influence contributions.
#' @param verbose Whether to print progress during fitting.
#' @details All penalties are zero by default. Regression ridge penalties apply to the
#' coefficients on the original covariate scale. Positive penalties change the
#' fitted working models. The fitted object records the chosen penalties and convergence
#' information. Weight truncation changes the numerical estimator and is therefore
#' an explicit option.
#'
#' The coefficient penalty is one half of the specified value times the sum of
#' squared coefficients, added to the negative summed likelihood or partial
#' likelihood. In a dependent joint copula fit, \code{association_ridge} penalizes
#' the squared distance of the copula parameter from its independence limit,
#' and \code{baseline_penalty} multiplies the sum of the final cumulative baseline
#' hazards over outcomes. These latter two penalties apply only to dependent
#' joint copula fits. Separate Cox margins under independence use coefficient
#' ridge only. Kaplan--Meier models have no regression coefficients to penalize.
#'
#' Let \eqn{m} be the number of independent randomization units. With
#' \code{finite_sample = TRUE} and \eqn{m > 2}, the empirical kernel projection
#' is multiplied by \eqn{(m-1)/(m-2)} before the contribution from fitting the
#' working models is added. Under independent randomization, the covariance
#' uses the centered sum of outer products divided by \eqn{m(m-1)}; with
#' \code{FALSE}, the denominator is \eqn{m^2}. Under complete randomization,
#' contributions are centered within arms. With the adjustment, the sum for
#' arm \eqn{a} is multiplied by \eqn{m_a/(m_a-1)}, where \eqn{m_a} is its
#' number of units, before the combined sum is divided by \eqn{m^2}. These
#' choices affect standard errors and intervals, while retaining the same
#' point estimates.
#' @return A validated control list of class \code{"drwins_control"}.
#' @examples drwins_control()
#' drwins_control(ridge = 0.001, maxit = 2000)
#' drwins_control(censor_ridge = 0.001, outcome_ridge = 0)
#' @seealso \code{\link{drwins}}
#' @export
drwins_control <- function(
    ridge = 0, censor_ridge = ridge, outcome_ridge = ridge,
    association_ridge = 0, baseline_penalty = 0,
    maxit = 1000, tol = 1e-7, finite_sample = TRUE,
    weight_floor = 1e-8, positivity = c("error", "truncate"),
    max_states = 200000, max_bytes = 5e8, keep_models = TRUE,
    keep_influence = TRUE, verbose = FALSE) {
  scalar_number <- function(x, name, lower = 0, strict = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || if (strict) x <= lower else x < lower) {
      stop(name, " must be a finite numeric scalar ",
           if (strict) "greater than " else "at least ", lower, ".",
           call. = FALSE)
    }
    unname(as.numeric(x))
  }
  scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    x
  }
  ridge <- scalar_number(ridge, "ridge")
  censor_ridge <- scalar_number(censor_ridge, "censor_ridge")
  outcome_ridge <- scalar_number(outcome_ridge, "outcome_ridge")
  association_ridge <- scalar_number(association_ridge, "association_ridge")
  baseline_penalty <- scalar_number(baseline_penalty, "baseline_penalty")
  maxit <- scalar_number(maxit, "maxit", 0, TRUE)
  if (maxit != floor(maxit) || maxit > .Machine$integer.max)
    stop("maxit must be a positive integer.", call. = FALSE)
  tol <- scalar_number(tol, "tol", 0, TRUE)
  weight_floor <- scalar_number(weight_floor, "weight_floor", 0, TRUE)
  if (weight_floor >= 1)
    stop("weight_floor must be less than one.", call. = FALSE)
  max_states <- scalar_number(max_states, "max_states", 0, TRUE)
  if (max_states != floor(max_states) || max_states > .Machine$integer.max)
    stop("max_states must be a positive whole number no larger than ",
         ".Machine$integer.max.", call. = FALSE)
  max_bytes <- scalar_number(max_bytes, "max_bytes", 0, TRUE)
  ans <- list(
    ridge = ridge, censor_ridge = censor_ridge,
    outcome_ridge = outcome_ridge, association_ridge = association_ridge,
    baseline_penalty = baseline_penalty, maxit = as.integer(maxit),
    tol = tol, finite_sample = scalar_flag(finite_sample, "finite_sample"),
    weight_floor = weight_floor, positivity = match.arg(positivity),
    max_states = max_states, max_bytes = max_bytes,
    keep_models = scalar_flag(keep_models, "keep_models"),
    keep_influence = scalar_flag(keep_influence, "keep_influence"),
    verbose = scalar_flag(verbose, "verbose")
  )
  class(ans) <- "drwins_control"
  ans
}
