# Fit the two-arm estimator. Native computation evaluates the bilinear
# comparison of augmented outcome measures; it does not sample pairs.
#' Fit Doubly Robust Win Statistics for Randomized Trials
#'
#' @name drwins
#' @rdname drwins
#' @description Estimate treatment win, loss, and tie probabilities at specified restriction
#' times for a prioritized composite of terminal and nonterminal events. The
#' function combines outcome regression with inverse probability of censoring
#' weighting and reports the net benefit, win ratio, and win odds. It accepts
#' longitudinal event records, fits the working models separately in the two
#' treatment arms, and provides analytic standard errors and confidence intervals.
#' The same interface handles individually randomized trials and
#' cluster-randomized trials.
#' @param data A nonempty data frame in longitudinal event format. Each row is an
#' observed nonterminal event or the participant's final observation. Each
#' participant must have exactly one final row with status 0 or 1. Required
#' columns and baseline covariates are described under \sQuote{Details} and in
#' \code{\link{drwins-data-format}}.
#' @param outcome A one-sided formula, such as \code{~ age + sex}, used for every
#' outcome component, or a list containing one one-sided formula per priority.
#' An unnamed list follows the order of \code{priorities}. A named list may use
#' all priority codes, such as \code{"1"} and \code{"2"}, or all their labels,
#' such as \code{"Death"} and \code{"Hospitalization"}. The default is
#' \code{~ 1}. Models use baseline covariates and are fitted separately by
#' arm, so the treatment column is excluded from the formula.
#' @param censoring A one-sided formula for the censoring model, fitted separately
#' in each treatment arm. The default is \code{~ 1}. The covariates can differ
#' from those in \code{outcome}.
#' @param id A single character string naming the participant identifier column.
#' The default is \code{"id"}. An identifier can appear on several event rows.
#' In a cluster-randomized trial, the combination of cluster and participant
#' identifiers identifies a participant, so local participant identifiers may
#' repeat across clusters.
#' @param time A single character string naming the time column. The default is
#' \code{"time"}. Values must be finite, nonnegative numbers measured from a
#' common time origin, using the same units as \code{times} and
#' \code{thresholds}.
#' @param status A single character string naming the event-code column. The
#' default is \code{"status"}. Codes are nonnegative integers: 0 indicates the
#' end of observation without the terminal event, 1 is the terminal event, and
#' 2 and higher identify nonterminal event types. A final status of 0 can
#' represent either dropout or administrative end of follow-up.
#' @param treatment A single character string naming the treatment column.
#' The default is \code{"treatment"}. Exactly two groups must be observed.
#' Assignment must be constant within a participant and within a cluster.
#' @param cluster A single character string naming the cluster identifier
#' column, or \code{NULL}, the default. With \code{NULL}, each participant is
#' an independent randomization unit. Supplying a cluster column uses clusters
#' as the independent units for estimation and variance calculation. A dataset
#' in which every cluster has one participant is handled as an IRT.
#' @param outcome_model The model for the component event times. The default,
#' \code{"cox"}, uses proportional hazards margins with baseline covariates
#' and the selected survival copula. \code{"nonparametric"} uses a joint model
#' with covariate-free margins and requires intercept-only outcome formulas.
#' \code{"km"} uses Kaplan--Meier working margins and requires intercept-only
#' formulas and the independence copula. See \sQuote{Working models} below.
#' The default \code{outcome = ~ 1} does not change \code{"cox"} to
#' \code{"km"}.
#' @param censoring_model \code{"cox"}, the default, for a Cox proportional
#' hazards censoring model, or \code{"km"} for a Kaplan--Meier censoring model.
#' The Kaplan--Meier option requires \code{censoring = ~ 1}.
#' @param copula A survival copula describing dependence among the component
#' event times of one participant. Choices are \code{"gumbel"}, the default,
#' \code{"independence"}, \code{"clayton"}, \code{"frank"}, and
#' \code{"joe"}. Association parameters for the four dependent families are
#' estimated separately by arm over their positive dependence ranges. With
#' \code{outcome_model = "km"}, omitting \code{copula} selects independence;
#' explicitly requesting a dependent family gives an error. A single-endpoint
#' analysis also uses independence automatically.
#' @param times A numeric vector of finite, strictly positive restriction
#' times. Explicit times are sorted and duplicates are removed. Every time must
#' be no greater than the maximum observed follow-up in either arm. If
#' \code{NULL}, the default, times are selected using \code{time_probs} from
#' participants' final observed follow-up. See \sQuote{Restriction times}.
#' @param time_probs A numeric vector of probabilities strictly between zero
#' and one used when \code{times = NULL}. The default \code{c(.25, .5, .75)}
#' selects the first quartile, median, and third quartile. Quantiles use one final
#' follow-up time per participant, including in CRTs. This argument must contain
#' valid probabilities even when explicit \code{times} are supplied.
#' @param target \code{"both"}, the default, \code{"cluster"}, or
#' \code{"individual"}. Cluster-average probabilities give equal weight to
#' clusters; individual-average probabilities give equal weight to participants.
#' With singleton clusters, \code{"both"} returns the individual-average
#' result once because the two targets coincide. An explicit
#' \code{target = "cluster"} is also allowed for singleton clusters.
#' @param allocation \code{"independent"}, the default, for independent
#' Bernoulli assignment with known treatment probabilities, or
#' \code{"complete"} for simple complete randomization with fixed numbers of
#' randomization units in each arm. The randomization unit is the participant
#' in an IRT and the cluster in a CRT.
#' @param treatment_probability The known probability of assignment to the
#' intervention group under independent randomization. The default is 0.5.
#' Supply a scalar for a common probability or one value per randomization unit.
#' All values must lie strictly between zero and one. Names, when supplied,
#' must identify every randomization unit exactly once. An unnamed vector
#' follows sorted character identifiers; a named vector is preferable when
#' probabilities differ. Under complete randomization, supply a valid scalar;
#' the calculation uses the fixed observed arm totals.
#' @param treatment_levels The two observed treatment labels in control,
#' intervention order. For example, \code{c("Usual care", "Intervention")}
#' makes the latter group the first member of every clinical comparison.
#' With \code{NULL}, factor levels retain their existing order after unused
#' levels are dropped, and other labels are sorted as character strings.
#' Set this argument explicitly to control the direction of comparison.
#' The label \code{"Overall"} is reserved for censoring summaries.
#' @param priorities A vector of unique positive event codes in clinical
#' priority order, with terminal code 1 first. Optional unique names provide
#' labels, for example \code{c(Death = 1L, Hospitalization = 2L)}.
#' With \code{NULL}, the function uses code 1 and every observed positive code
#' in increasing order. Every positive code present in \code{data} must appear
#' in this vector. A code can be supplied even if no corresponding event is
#' observed. At most 20 priorities are supported, subject to the computational
#' limits in \code{control}.
#' @param thresholds Nonnegative indifference thresholds in the same time
#' units as the outcomes. Supply one value for all priorities or one value per
#' priority in its specified order. The default is zero. A difference no greater
#' than the threshold on a priority leads to comparison on the next priority.
#' @param recurrent \code{"first"}, the default, uses the first occurrence of
#' each nonterminal event type and reports a warning when later records are
#' removed. \code{"error"} rejects repeated records of the same nonterminal
#' event type. Both options concern first-event times, rather than a comparison
#' based on recurrent-event counts.
#' @param conf_level A number strictly between zero and one specifying the
#' confidence level stored in the fit. The default is 0.95. After fitting,
#' \code{summary(fit, conf.level = ...)} and \code{plot(fit, level = ...)}
#' can use another confidence level without refitting the models.
#' @param inference \code{"sandwich"}, the default, calculates analytic
#' standard errors including estimation of the working models and, for the
#' individual-average target in CRTs, mean cluster size. \code{"none"}
#' calculates point estimates only. Sandwich inference requires at least three
#' independent units; complete randomization additionally requires at least two
#' units in each arm.
#' @param control A control object returned by \code{\link{drwins_control}}, or
#' a named list of arguments to that function. It specifies optional penalties,
#' optimization settings, variance adjustments, numerical limits, and storage.
#' All penalties are zero by default. For example,
#' \code{drwins_control(ridge = 0.001)} explicitly requests coefficient ridge
#' penalties for both the outcome and censoring models.
#' @details \subsection{Longitudinal event records}{
#' For a composite of death and hospitalization, a hospitalization at time 0.7
#' followed by death at time 1.8 is represented by a status-2 row at 0.7 and a
#' status-1 row at 1.8. Hospitalization followed by the end of observation at
#' 1.8 instead has a status-0 final row. A participant without an observed
#' nonterminal event has only the final row. Row order does not affect fitting.
#'
#' Each participant must have exactly one final status-0 or status-1 row.
#' Events recorded after that final time are rejected. A nonterminal event
#' recorded at the terminal-event time is omitted with a warning, since it
#' does not precede death. Baseline covariates in the formulas must be complete
#' and constant across the participant's rows. Missing covariates cause an
#' error; the function does not remove participants automatically.
#'
#' Formulas can contain factors, interactions, and transformations such as
#' \code{I(x^2)} or \code{log(cluster_size)}. The columns designated as
#' \code{id}, \code{time}, \code{status}, \code{treatment}, or \code{cluster}
#' are excluded from the formulas. Offsets and time-varying covariates are not
#' supported. Cluster-level covariates and a separately recorded baseline
#' cluster size can be included in a CRT formula. Cluster size is not added
#' automatically as a covariate.
#' }
#' \subsection{Clinical comparison and win statistics}{
#' At a restriction time \eqn{\tau}, component event times are restricted to
#' \eqn{\tau}. An event that has not occurred by that time receives value
#' \eqn{\tau}; the same convention applies to a nonterminal event precluded
#' by death. The first participant in a comparison is from the intervention
#' group, as determined by \code{treatment_levels}. Longer times are preferred.
#' The terminal event is compared first. A difference exceeding its threshold
#' determines the result; otherwise the comparison proceeds through the
#' remaining priorities. A tie occurs when all differences are within their
#' thresholds. With zero thresholds, an exact tie on the terminal event still
#' allows a lower-priority event to determine the comparison.
#'
#' Let \eqn{p_+}{p_win}, \eqn{p_-}{p_loss}, and \eqn{p_=}{p_tie} denote the
#' win, loss, and tie probabilities. The reported net benefit is
#' \eqn{\mathrm{NB} = p_+ - p_-}{NB = p_win - p_loss}; the win ratio is
#' \eqn{\mathrm{WR} = p_+/p_-}{WR = p_win / p_loss}; and the win odds are
#' \eqn{\mathrm{WO} = (p_+ + p_=/2)/(p_- + p_=/2)}{WO = (p_win + p_tie/2) / (p_loss + p_tie/2)}.
#' Equivalently, \eqn{\mathrm{WO} = (1+\mathrm{NB})/(1-\mathrm{NB})}{WO = (1 + NB) / (1 - NB)}.
#' A positive net benefit, win ratio above one, or win odds above one favors
#' the intervention. The tie probability is estimated as one minus the sum
#' of the estimated win and loss probabilities.
#' }
#' \subsection{Target population and randomization}{
#' An IRT compares two independent participants, one evaluated under each
#' condition. In a CRT, comparisons involve participants from different
#' clusters. The cluster-average target first selects clusters with equal
#' weight and then selects a participant within each cluster with equal
#' weight. The individual-average target selects participants with equal
#' weight, giving larger clusters more weight. The two targets can differ
#' when cluster size is associated with outcomes or treatment effects.
#'
#' The independent randomization option uses the supplied known assignment
#' probabilities. Complete randomization uses the numbers of independent
#' units in each arm and centers the full influence contributions within
#' arms for variance estimation. These options implement independent and
#' simple complete randomization; stratified or matched designs require
#' variance calculations appropriate to those designs.
#' }
#' \subsection{Working models}{
#' Outcome and censoring models are fitted separately within each treatment
#' arm. With dependent event times, the Cox outcome option joins proportional
#' hazards margins through the selected survival copula and estimates the
#' regression coefficients, baseline hazards, and association parameter using
#' the observed event histories. The \code{"nonparametric"} option fits a
#' covariate-free joint model with the selected copula. Under independence,
#' the component models are fitted separately.
#'
#' The Kaplan--Meier outcome option uses product-limit estimates and an
#' independence model for the event times. Death acts as censoring for a
#' nonterminal-event margin in that option. A dependent, covariate-free
#' outcome model is specified instead by \code{outcome_model = "nonparametric"}
#' and a dependent copula. The censoring model can be Cox or Kaplan--Meier
#' regardless of the chosen outcome model, subject to their formula requirements.
#'
#' The copula describes dependence among outcomes of the same participant.
#' Within-cluster dependence is accounted for in the cluster-level covariance
#' calculation. The implemented marginal working models can include baseline
#' information from the cluster; they do not estimate a joint distribution of
#' all participants' event histories within a cluster.
#'
#' The doubly robust estimator is consistent if the censoring models are
#' correctly specified in both arms or the joint outcome models are correctly
#' specified in both arms, under the sampling, randomization, covariate dependent
#' censoring, positivity, and regularity assumptions. The outcome-model condition
#' concerns both the margins and their dependence. Conditional independence of
#' censoring and the event outcomes given treatment and the supplied baseline
#' information remains an identifying assumption.
#' }
#' \subsection{Restriction times and censoring rates}{
#' Explicit restriction times define the follow-up windows for the clinical
#' comparisons. With \code{times = NULL}, the function uses type-7 sample
#' quantiles of participants' final observed follow-up times. Each participant
#' contributes once, regardless of the number of event rows or cluster size.
#' Duplicate or nonpositive default quantiles are removed with a warning, so
#' the result can contain fewer times than requested. The actual fitted times
#' are stored in \code{fit$times}. Default quantiles provide exploratory
#' choices; prespecified analyses should supply their clinical restriction
#' times explicitly. Models are fitted using observations through
#' \code{max(times)}.
#'
#' Censoring before \eqn{\tau} is counted when a participant's final status
#' is 0 and the final time is strictly less than \eqn{\tau}. A final status-0
#' row at \eqn{\tau} or later provides observation through that time. Death
#' before \eqn{\tau} also resolves the restricted outcomes. Consequently,
#' the fraction of all final rows with status 0 is not the censoring rate
#' at a particular restriction time.
#'
#' The censoring table reports each treatment arm and the overall sample
#' at every target and fitted time. For the individual-average target,
#' \code{rate} is the participant proportion censored. For the cluster-average
#' target it is the mean of the participant censoring proportions within
#' clusters. Counts in \code{n} and \code{censored} refer to participants;
#' the cluster-average \code{rate} need not equal \code{censored / n}.
#' }
#' \subsection{Inference and numerical controls}{
#' The sandwich covariance accounts for reuse of participants or clusters
#' across pairwise comparisons and for fitting the working models. It uses
#' participants as independent units in IRTs and clusters in CRTs. For the
#' individual-average target, it includes estimation of mean cluster size
#' and its covariance with the estimated comparison probabilities.
#'
#' Probability and net benefit intervals use their natural scales. Win ratio
#' and win odds intervals are calculated on the logarithmic scale and
#' exponentiated. Intervals are pointwise at each fitted restriction time.
#' The reported two-sided Wald tests use null value zero for net benefit and
#' one for the two ratios. Changing the confidence level in a summary or plot
#' recalculates intervals from the stored standard errors.
#'
#' Raw augmented estimates are retained, including estimates outside the
#' probability bounds in a finite sample. Undefined ratios or logarithmic
#' intervals are returned as missing. Numerical convergence, censoring weights,
#' and probability-bound diagnostics are available in \code{fit$diagnostics}.
#' Ridge penalties and weight truncation are explicit choices in
#' \code{\link{drwins_control}}. The exact outcome-support calculation retains
#' the specified priorities; support or memory requirements exceeding its
#' configured limits produce an error.
#' }
#' @return A list of class \code{"drwins"}, with the following components.
#' \describe{
#' \item{call}{The original function call.}
#' \item{estimates}{A data frame with one row per target, time, and measure.
#' Its columns are described in \sQuote{Estimate table} below.}
#' \item{covariance}{The joint covariance matrix of the estimated win and loss
#' probabilities across all fitted targets and times. Its row and column names
#' use \code{target:time:probability}. The matrix includes the sample-size
#' scaling for the estimates. Under \code{inference = "none"}, its entries
#' are missing.}
#' \item{influence}{A matrix of estimated influence contributions, with one
#' row per independent participant or cluster and columns matching
#' \code{covariance}. It is \code{NULL} if \code{control$keep_influence = FALSE}
#' or \code{inference = "none"}.}
#' \item{nuisance_influence}{The part of the influence contribution due to
#' fitting the working models, with the same rows and columns as
#' \code{influence}. Subject to the same storage controls.}
#' \item{censoring}{A data frame with columns \code{target}, \code{time},
#' \code{arm}, \code{n}, \code{clusters}, \code{censored}, and \code{rate}.
#' The rate is stored as a fraction; summaries print it as a percentage.}
#' \item{design}{Trial information, including \code{type} (\code{"IRT"}
#' or \code{"CRT"}), the participant and independent-unit counts \code{n}
#' and \code{m}, counts by arm, cluster sizes and labels, treatment levels,
#' column mappings, allocation, and information about the selected times.}
#' \item{models}{The fitted arm-specific outcome and censoring model lists,
#' including coefficients, baseline survival estimates, and outcome dependence
#' parameters where applicable. This component is \code{NULL} when
#' \code{control$keep_models = FALSE}.}
#' \item{diagnostics}{A list containing data-processing information, model
#' convergence information, censoring-weight summaries, outcome-support
#' summaries, and indicators of estimates outside the probability bounds.}
#' \item{settings}{The model choices, normalized formulas, priorities and
#' labels, thresholds, treatment levels, confidence level, inference option,
#' allocation, and numerical controls used for the fit.}
#' \item{times}{The sorted, distinct restriction times actually fitted.}
#' \item{target}{The target names represented in the result.}
#' \item{formula}{A list containing the normalized outcome formulas and
#' censoring formula.}
#' \item{version}{The package version used to create the object.}
#' }
#' Use \code{\link{summary.drwins}} to select and display results,
#' \code{\link{plot.drwins}} to plot them, or the extraction methods documented
#' in \code{\link{drwins-methods}}.
#' @section Estimate table:
#' \code{fit$estimates} contains the following columns.
#' \describe{
#' \item{target, time}{The target population and restriction time.}
#' \item{measure}{\code{"win"}, \code{"loss"}, \code{"tie"},
#' \code{"NB"}, \code{"WR"}, or \code{"WO"}.}
#' \item{estimate}{The point estimate on its natural scale, including the
#' natural ratio scale for WR and WO.}
#' \item{se}{The standard error on the natural scale.}
#' \item{se_log}{The standard error of the logarithm for WR and WO;
#' missing for the other measures.}
#' \item{lower, upper}{Confidence limits on the natural scale. Ratio limits
#' have already been exponentiated.}
#' \item{p_value}{A two-sided normal Wald p-value for NB, WR, or WO. Missing
#' for the three probabilities or when the relevant test is undefined.}
#' \item{scale}{\code{"natural"} for probabilities and NB, or
#' \code{"log"} for WR and WO, identifying the scale used to construct
#' intervals. This column does not mean that \code{estimate} or \code{se}
#' has been log-transformed.}
#' }
#' With \code{inference = "none"}, standard errors, confidence limits, and
#' p-values are missing; the point estimates remain available.
#' @seealso \code{\link{summary.drwins}}, \code{\link{plot.drwins}},
#' \code{\link{drwins-methods}}, \code{\link{drwins_control}},
#' \code{\link{drwins-data-format}}, \code{\link{example_irt}},
#' \code{\link{example_crt}}. The introductory vignette is available with
#' \code{vignette("DRwinStats")}.
#' @examples ## 1. Individually randomized trial with default restriction times.
#' data(example_irt)
#' final_irt <- example_irt[example_irt$status %in% c(0L, 1L), ]
#' ids <- unlist(lapply(split(final_irt$id, final_irt$treatment),
#'                     head, 20L), use.names = FALSE)
#' small_irt <- example_irt[example_irt$id %in% ids, ]
#'
#' # Keep all event rows for 20 participants from each arm.
#' fit_irt <- drwins(
#'   small_irt, outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence",
#'   priorities = c(Death = 1L, Hospitalization = 2L),
#'   allocation = "complete"
#' )
#' fit_irt$times                       # Default 25th, 50th, 75th percentiles.
#' summary(fit_irt, measure = c("WR", "NB"))
#' fit_irt$censoring                   # Rates are stored as fractions.
#' plot(fit_irt, measure = "WR")
#'
#' ## 2. Different outcome formulas, Cox censoring, and a Gumbel copula.
#' ids_model <- unlist(lapply(split(final_irt$id, final_irt$treatment),
#'                           head, 40L), use.names = FALSE)
#' model_irt <- example_irt[example_irt$id %in% ids_model, ]
#' fit_model <- drwins(
#'   model_irt,
#'   outcome = list(Death = ~ age + x, Hospitalization = ~ age),
#'   censoring = ~ age + x,
#'   outcome_model = "cox", censoring_model = "cox",
#'   copula = "gumbel",
#'   priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 2), allocation = "complete"
#' )
#' summary(fit_model, measure = c("WR", "WO", "NB"), times = 2)
#' coef(fit_model, measure = "WR")
#' confint(fit_model, measure = "WR", level = 0.90)
#'
#' ## 3. Cluster-randomized trial, retaining whole clusters.
#' data(example_crt)
#' cluster_arm <- unique(example_crt[c("cluster", "treatment")])
#' keep_clusters <- unlist(
#'   lapply(split(cluster_arm$cluster, cluster_arm$treatment), head, 4L),
#'   use.names = FALSE
#' )
#' small_crt <- example_crt[example_crt$cluster %in% keep_clusters, ]
#' fit_crt <- drwins(
#'   small_crt, cluster = "cluster",
#'   outcome = ~ 1, censoring = ~ 1,
#'   outcome_model = "km", censoring_model = "km",
#'   copula = "independence",
#'   priorities = c(Death = 1L, Hospitalization = 2L),
#'   times = c(1, 2), target = "both", allocation = "complete"
#' )
#' summary(fit_crt, measure = "WR")
#' summary(fit_crt, measure = "NB", target = "individual", times = 2)
#' plot(fit_crt, measure = "WR")
#' head(as.data.frame(fit_crt))
#'
#' # Optional controls can be passed as control = drwins_control(...).
#' drwins_control(censor_ridge = 0.001, outcome_ridge = 0)
#' @author Xi Fang \email{xfang@mcw.edu}
#' @keywords models
#' @concept win statistics
#' @concept cluster-randomized trials
#' @export
drwins <- function(
    data, outcome = ~ 1, censoring = ~ 1,
    id = "id", time = "time", status = "status", treatment = "treatment",
    cluster = NULL, outcome_model = c("cox", "km", "nonparametric"),
    censoring_model = c("cox", "km"),
    copula = c("gumbel", "independence", "clayton", "frank", "joe"),
    times = NULL, time_probs = c(.25, .5, .75),
    target = c("both", "cluster", "individual"),
    allocation = c("independent", "complete"), treatment_probability = .5,
    treatment_levels = NULL, priorities = NULL, thresholds = 0,
    recurrent = c("first", "error"), conf_level = .95,
    inference = c("sandwich", "none"), control = drwins_control()) {
  call <- match.call()
  copula_supplied <- !missing(copula)
  outcome_model <- match.arg(outcome_model)
  censoring_model <- match.arg(censoring_model)
  copula <- match.arg(copula)
  target <- match.arg(target)
  allocation <- match.arg(allocation)
  recurrent <- match.arg(recurrent)
  inference <- match.arg(inference)
  conf_level <- .drw_conf_level(conf_level)
  if (!is.list(control)) stop("control must be a drwins_control() object or a named list.", call. = FALSE)
  if (!inherits(control, "drwins_control")) control <- do.call(drwins_control, control)
  prepared <- .drw_prepare_data(
    data = data, outcome = outcome, censoring = censoring, id = id,
    time = time, status = status, treatment = treatment, cluster = cluster,
    outcome_model = outcome_model, censoring_model = censoring_model,
    times = times, time_probs = time_probs, target = target,
    allocation = allocation, treatment_probability = treatment_probability,
    treatment_levels = treatment_levels, priorities = priorities,
    thresholds = thresholds, recurrent = recurrent)
  prepared$horizon <- max(prepared$times)
  if (ncol(prepared$event_times) == 1L || (outcome_model == "km" && !copula_supplied))
    copula <- "independence"
  if (outcome_model == "km" && copula != "independence")
    stop("Kaplan-Meier outcome margins require copula = 'independence'. ",
         "For a copula with covariate-free joint margins, use ",
         "outcome_model = 'nonparametric' and outcome = ~ 1.", call. = FALSE)
  support_end <- vapply(0:1, function(a)
    max(prepared$stop[prepared$arm == a]), numeric(1))
  if (any(prepared$times > min(support_end)))
    stop("A requested restriction time exceeds the maximum observed follow-up ",
         "in at least one treatment arm. Choose supported times explicitly.", call. = FALSE)
  if (inference == "sandwich") {
    if (prepared$m < 3L)
      stop("Analytic variance requires at least three independent units. ",
           "Use inference = 'none' for point estimates only.", call. = FALSE)
    if (allocation == "complete" && any(table(prepared$cluster_arm) < 2L))
      stop("Within-arm variance estimation requires at least two independent ",
           "units in each arm. Use inference = 'none' for point estimates only.", call. = FALSE)
  }
  if (control$verbose) message("Fitting outcome and censoring models in each treatment arm.")
  models <- .drw_fit_models(prepared, outcome_model, censoring_model, copula,
                           control, inference)
  models$compute_if <- inference == "sandwich"
  if (control$verbose) message("Evaluating win probabilities at ",
                               paste(signif(prepared$times, 5), collapse = ", "), ".")
  target_code <- as.integer(prepared$targets == "individual")
  result <- .Call("drw_estimate", prepared, models, target_code, control,
                  PACKAGE = "DRwinStats")
  m <- prepared$m
  index <- do.call(rbind, lapply(prepared$times, function(t)
    data.frame(time = t, target = prepared$targets, stringsAsFactors = FALSE)))
  probability <- result$estimate
  influence <- result$projection
  nuisance <- matrix(0, m, ncol(influence))
  if (inference == "sandwich") {
    nuisance_if <- .drw_model_influence(models, m)
    if (ncol(nuisance_if) != ncol(result$derivative))
      stop("Internal nuisance derivative dimensions are inconsistent.", call. = FALSE)
    nuisance <- nuisance_if %*% t(result$derivative)
    influence <- influence + nuisance
  }
  mean_size <- mean(prepared$cluster_size)
  for (j in seq_len(nrow(index))) if (index$target[j] == "individual") {
    cols <- 2L * j + c(-1L, 0L)
    probability[j, ] <- probability[j, ] / mean_size^2
    influence[, cols] <- influence[, cols, drop = FALSE] / mean_size^2 -
      outer(prepared$cluster_size - mean_size, 2 * probability[j, ] / mean_size)
    nuisance[, cols] <- nuisance[, cols, drop = FALSE] / mean_size^2
  }
  labels <- unlist(lapply(seq_len(nrow(index)), function(j)
    paste(index$target[j], format(index$time[j], digits = 12, trim = TRUE),
          c("win", "loss"), sep = ":")), use.names = FALSE)
  colnames(influence) <- colnames(nuisance) <- labels
  rownames(influence) <- rownames(nuisance) <- prepared$design$cluster_labels
  if (inference == "sandwich") {
    covariance <- .drw_covariance(influence, prepared$cluster_arm,
                                  allocation, control$finite_sample)
  } else covariance <- matrix(NA_real_, ncol(influence), ncol(influence),
                               dimnames = list(labels, labels))
  estimates <- .drw_summaries(probability, covariance, index, conf_level)
  outside <- probability[, 1L] < -1e-8 | probability[, 2L] < -1e-8 |
    rowSums(probability) > 1 + 1e-8
  if (any(outside)) warning(
    "Augmented probability estimates are outside their probability bounds ",
    "for ", sum(outside), " target-time combination(s). Raw estimates are ",
    "retained; ratio intervals are returned only when defined.", call. = FALSE)
  if (any(result$truncated_weights > 0)) warning(
    "Censoring weights were truncated as explicitly requested. ",
    "See diagnostics$weights for the affected times.", call. = FALSE)
  diagnostics <- list(
    data = prepared$diagnostics, models = models$diagnostics,
    weights = data.frame(time = prepared$times,
                         min_survival = result$min_censor_survival,
                         max_weight = result$max_censor_weight,
                         truncated = result$truncated_weights),
    support = data.frame(time = prepared$times, states = result$states,
                         augmentation_mass_error = result$mass_error),
    outside_probability_bounds = cbind(index, outside = outside))
  fit <- list(
    call = call, estimates = estimates, covariance = covariance,
    influence = if (control$keep_influence && inference == "sandwich") influence else NULL,
    nuisance_influence = if (control$keep_influence && inference == "sandwich") nuisance else NULL,
    censoring = prepared$censoring, design = prepared$design,
    models = if (control$keep_models) models else NULL,
    diagnostics = diagnostics,
    settings = list(outcome_model = outcome_model, censoring_model = censoring_model,
                    copula = copula, outcome_formula = prepared$outcome_formulas,
                    censoring_formula = prepared$censoring_formula,
                    outcome = prepared$outcome_formulas, censoring = prepared$censoring_formula,
                    priorities = prepared$priorities, labels = prepared$labels,
                    thresholds = prepared$thresholds, conf_level = conf_level,
                    inference = inference, control = control,
                    treatment_levels = prepared$treatment_levels,
                    default_times = prepared$design$default_times,
                    time_probs = time_probs, allocation = allocation),
    times = prepared$times, target = prepared$targets,
    formula = list(outcome = prepared$outcome_formulas, censoring = prepared$censoring_formula),
    version = "0.1.2")
  class(fit) <- "drwins"
  fit
}

.drw_model_influence <- function(models, m) {
  values <- list()
  put <- function(x, p) {
    if (is.null(x)) x <- matrix(0, m, p)
    x <- as.matrix(x)
    if (!identical(dim(x), c(as.integer(m), as.integer(p))))
      stop("An internal nuisance influence matrix has incorrect dimensions.", call. = FALSE)
    if (any(!is.finite(x)))
      stop("A nuisance influence matrix contains nonfinite values.", call. = FALSE)
    values[[length(values) + 1L]] <<- x
  }
  for (a in 1:2) {
    z <- models$outcome[[a]]
    for (j in seq_along(z$beta)) {
      put(z$beta_if[[j]], length(z$beta[[j]]))
      put(z$survival_if[[j]], length(z$base_survival[[j]]))
    }
    put(matrix(z$theta_if, nrow = m), 1L)
  }
  for (a in 1:2) {
    z <- models$censoring[[a]]
    put(z$beta_if, length(z$beta))
    put(z$survival_if, length(z$base_survival))
  }
  do.call(cbind, values)
}

.drw_covariance <- function(influence, arm, allocation, correction) {
  m <- nrow(influence)
  V <- matrix(0, ncol(influence), ncol(influence),
              dimnames = list(colnames(influence), colnames(influence)))
  if (allocation == "complete") {
    for (a in 0:1) {
      z <- influence[arm == a, , drop = FALSE]
      z <- sweep(z, 2L, colMeans(z))
      V <- V + crossprod(z) * (if (correction) nrow(z)/(nrow(z) - 1L) else 1)
    }
    V <- V / m^2
  } else {
    z <- sweep(influence, 2L, colMeans(influence))
    V <- crossprod(z) / (if (correction) m * (m - 1L) else m^2)
  }
  V
}

.drw_summaries <- function(probability, covariance, index, conf_level) {
  z <- stats::qnorm((1 + conf_level)/2)
  out <- vector("list", nrow(index))
  for (i in seq_len(nrow(index))) {
    p <- probability[i, 1L]; q <- probability[i, 2L]; nb <- p - q
    tie <- 1 - p - q
    wr <- if (q > 0 && p >= 0) p/q else NA_real_
    wo <- if (abs(nb) < 1) (1 + nb)/(1 - nb) else NA_real_
    est <- c(p, q, tie, nb, wr, wo)
    g <- rbind(c(1, 0), c(0, 1), c(-1, -1), c(1, -1),
               if (is.finite(wr)) c(1/q, -p/q^2) else c(NA, NA),
               if (is.finite(wo)) 2*c(1, -1)/(1 - nb)^2 else c(NA, NA))
    cols <- 2L*i + c(-1L, 0L)
    V <- covariance[cols, cols, drop = FALSE]
    vv <- rowSums((g %*% V) * g)
    se <- sqrt(pmax(vv, 0))
    se_log <- rep(NA_real_, 6L)
    if (p > 0 && q > 0) {
      gl <- c(1/p, -1/q)
      se_log[5L] <- sqrt(max(0, as.numeric(t(gl) %*% V %*% gl)))
    }
    if (abs(nb) < 1) {
      gl <- 2*c(1, -1)/(1 - nb^2)
      se_log[6L] <- sqrt(max(0, as.numeric(t(gl) %*% V %*% gl)))
    }
    lower <- est - z*se; upper <- est + z*se
    lower[5:6] <- exp(log(est[5:6]) - z*se_log[5:6])
    upper[5:6] <- exp(log(est[5:6]) + z*se_log[5:6])
    pvalue <- rep(NA_real_, 6L)
    if (is.finite(se[4L]) && se[4L] > 0) pvalue[4L] <- 2*stats::pnorm(-abs(nb/se[4L]))
    for (j in 5:6) if (is.finite(se_log[j]) && se_log[j] > 0 && is.finite(est[j]) && est[j] > 0)
      pvalue[j] <- 2*stats::pnorm(-abs(log(est[j])/se_log[j]))
    out[[i]] <- data.frame(target = index$target[i], time = index$time[i],
                           measure = c("win", "loss", "tie", "NB", "WR", "WO"),
                           estimate = est, se = se, se_log = se_log,
                           lower = lower, upper = upper, p_value = pvalue,
                           scale = c(rep("natural", 4), rep("log", 2)),
                           stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
