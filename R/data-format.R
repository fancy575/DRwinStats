#' Longitudinal Event Data for DRwinStats
#'
#' @name drwins-data-format
#' @rdname drwins-data-format
#' @description The input to \code{\link{drwins}} is a data frame with one row for each
#' observed event or final end of follow-up.
#' @details The columns specified by \code{id}, \code{time}, \code{status}, and
#' \code{treatment} identify the participant, time since a common origin, event
#' code, and randomized treatment. For a CRT, \code{cluster} additionally names
#' the cluster identifier. Participant IDs may repeat in different clusters.
#'
#' Code zero ends observation without the terminal event. Code one is the
#' terminal event and must be the highest-priority endpoint. Codes two and above
#' identify the lower-priority nonterminal endpoints. The \code{priorities}
#' argument determines their clinical order, with code one first; names supplied
#' on that vector become endpoint labels. At most 20 priorities are supported.
#' The state and memory limits can be reached with fewer priorities.
#'
#' Each participant requires exactly one final record with code zero or one.
#' For example, a hospitalization at time 0.5 followed by censoring at time 2
#' requires a code-two row at 0.5 and a code-zero row at 2. A hospitalization row
#' alone does not specify subsequent follow-up. A death at time 1.5 is instead
#' recorded with a code-one row at 1.5. There can be no records after that row.
#'
#' The analysis uses the first event of each nonterminal type. With
#' \code{recurrent = "first"}, later records of the same type are removed with
#' a warning. With \code{recurrent = "error"}, those records cause an error.
#' Nonterminal events recorded at death are removed with a warning, following
#' the requirement that a nonterminal event precede death.
#'
#' Time values must be finite and nonnegative. Treatment and every baseline
#' covariate used in a formula must be constant within each participant.
#' Treatment must also be constant within each cluster. Missing formula
#' covariates cause an error; participants are not deleted automatically.
#' Formulas are one-sided and can include factors, interactions, and ordinary
#' transformations. Models are fitted separately in the two treatment arms.
#'
#' The default restriction times are the 25th, 50th, and 75th percentiles
#' of the participant-level final observed follow-up times, using the default
#' R quantile definition (\code{type = 7}). Duplicate or nonpositive
#' quantiles are removed, with a warning. Supply \code{times} explicitly for
#' prespecified clinical restriction times.
#'
#' A final code-zero record can represent administrative end of follow-up.
#' At a requested restriction time \eqn{\tau}, the reported censoring rate counts
#' only participants whose code-zero record occurs strictly before \eqn{\tau}.
#' Participants with a terminal event before \eqn{\tau}, or who remain under
#' observation through \eqn{\tau}, have resolved restricted outcomes.
#' @examples long <- data.frame(
#'   id = c(1, 1, 2, 2, 3, 4),
#'   treatment = c(0, 0, 1, 1, 0, 1),
#'   time = c(0.5, 2, 0.8, 1.5, 2, 2),
#'   status = c(2, 0, 2, 1, 0, 0),
#'   age = c(65, 65, 71, 71, 59, 68)
#' )
#' long
#' @seealso \code{\link{drwins}}, \code{\link{example_irt}},
#'   \code{\link{example_crt}}
NULL
