.drw_column <- function(data, name, role, optional = FALSE) {
  if (optional && is.null(name)) return(NULL)
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(name) || !name %in% names(data))
    stop(role, " must name one column in data.", call. = FALSE)
  name
}

.drw_formula <- function(x, name) {
  if (!inherits(x, "formula") || length(x) != 2L)
    stop(name, " must be a one-sided formula, such as ~ age + sex.",
         call. = FALSE)
  x
}

.drw_outcome_formulas <- function(outcome, priorities, labels) {
  d <- length(priorities)
  if (d > 20L)
    stop("At most 20 outcome priorities are supported. The state and memory ",
         "limits can be reached with fewer priorities.", call. = FALSE)
  if (inherits(outcome, "formula")) {
    out <- rep(list(.drw_formula(outcome, "outcome")), d)
  } else {
    if (!is.list(outcome) || length(outcome) != d)
      stop("outcome must be a one-sided formula or one formula for each priority.",
           call. = FALSE)
    out <- outcome
    if (!is.null(names(out))) {
      if (anyNA(names(out)) || any(!nzchar(names(out))) ||
          anyDuplicated(names(out)))
        stop("Names of outcome formulas must be unique and nonempty.",
             call. = FALSE)
      wanted <- if (all(as.character(priorities) %in% names(out)))
        as.character(priorities) else labels
      if (!all(wanted %in% names(out)))
        stop("Name outcome formulas by the priority codes or their labels.",
             call. = FALSE)
      out <- out[wanted]
    }
    for (j in seq_len(d))
      out[[j]] <- .drw_formula(out[[j]], paste0("outcome[[", j, "]]"))
  }
  names(out) <- as.character(priorities)
  out
}

.drw_model_matrix <- function(formula, data, name) {
  frame <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.fail,
                       drop.unused.levels = TRUE),
    error = function(e)
      stop("Cannot construct the ", name, " model: ", conditionMessage(e),
           call. = FALSE)
  )
  tt <- attr(frame, "terms")
  if (length(attr(tt, "offset")))
    stop("Offsets are not supported in the ", name,
         " formula; include ordinary baseline covariates.", call. = FALSE)
  mm <- tryCatch(
    stats::model.matrix(tt, frame),
    error = function(e)
      stop("Cannot construct the ", name, " model matrix: ",
           conditionMessage(e), call. = FALSE)
  )
  if (nrow(mm) != nrow(data) || anyNA(mm) || any(!is.finite(mm)))
    stop("The ", name,
         " formula must give finite values for every participant.",
         call. = FALSE)
  intercept <- colnames(mm) == "(Intercept)"
  constant <- if (ncol(mm)) vapply(seq_len(ncol(mm)), function(j)
    length(unique(mm[, j])) <= 1L, logical(1)) else logical()
  omitted <- colnames(mm)[intercept | constant]
  mm <- mm[, !(intercept | constant), drop = FALSE]
  storage.mode(mm) <- "double"
  attr(mm, "dropped_columns") <- omitted
  attr(mm, "terms") <- tt
  mm
}

.drw_censoring_table <- function(stop, terminal, arm, cluster, times,
                                  targets, treatment_levels) {
  rows <- list()
  row <- 0L
  for (target in targets) for (tau in times) for (a in c(-1L, 0L, 1L)) {
    use <- if (a == -1L) rep(TRUE, length(stop)) else arm == a
    cens <- terminal[use] == 0L & stop[use] < tau
    cls <- cluster[use]
    rates <- if (target == "cluster")
      vapply(split(cens, cls), mean, numeric(1)) else cens
    row <- row + 1L
    rows[[row]] <- data.frame(
      target = target, time = tau,
      arm = if (a == -1L) "Overall" else treatment_levels[a + 1L],
      n = sum(use), clusters = length(unique(cls)),
      censored = sum(cens), rate = mean(rates),
      stringsAsFactors = FALSE
    )
  }
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

.drw_prepare_data <- function(
    data, outcome = ~ 1, censoring = ~ 1, id = "id", time = "time",
    status = "status", treatment = "treatment", cluster = NULL,
    outcome_model = "cox", censoring_model = "cox",
    times = NULL, time_probs = c(.25, .5, .75),
    target = "both", allocation = "independent",
    treatment_probability = .5, treatment_levels = NULL,
    priorities = NULL, thresholds = 0, recurrent = "first") {
  if (!is.data.frame(data) || !nrow(data))
    stop("data must be a nonempty data frame.", call. = FALSE)
  if (anyDuplicated(names(data)) || anyNA(names(data)))
    stop("Column names in data must be unique and nonmissing.", call. = FALSE)
  id <- .drw_column(data, id, "id")
  time <- .drw_column(data, time, "time")
  status <- .drw_column(data, status, "status")
  treatment <- .drw_column(data, treatment, "treatment")
  cluster <- .drw_column(data, cluster, "cluster", optional = TRUE)
  structural <- c(id, time, status, treatment, cluster)
  if (anyDuplicated(c(id, time, status, treatment)) ||
      (!is.null(cluster) && cluster %in% c(time, status, treatment)))
    stop("id, time, status, and treatment must name distinct columns; ",
         "cluster may coincide with id, but not with time, status, or treatment.",
         call. = FALSE)
  outcome_model <- match.arg(outcome_model, c("cox", "km", "nonparametric"))
  censoring_model <- match.arg(censoring_model, c("cox", "km"))
  allocation <- match.arg(allocation, c("independent", "complete"))
  target <- match.arg(target, c("both", "cluster", "individual"))
  recurrent <- match.arg(recurrent, c("first", "error"))
  tt <- data[[time]]
  if (!is.numeric(tt) || anyNA(tt) || any(!is.finite(tt)) || any(tt < 0))
    stop("Event and follow-up times must be finite, nonnegative numbers.",
         call. = FALSE)
  ee <- data[[status]]
  if (is.factor(ee)) ee <- as.character(ee)
  if (is.character(ee) && all(!is.na(ee) & grepl("^[0-9]+$", ee)))
    ee <- as.numeric(ee)
  if (!is.numeric(ee) || anyNA(ee) || any(!is.finite(ee)) ||
      any(ee < 0 | ee != floor(ee)))
    stop("status must contain nonnegative integer event codes.", call. = FALSE)
  if (any(ee > .Machine$integer.max))
    stop("Event codes exceed the integer range.", call. = FALSE)
  ee <- as.integer(ee)
  labels <- names(priorities)
  if (is.null(priorities)) priorities <- sort(unique(c(1L, ee[ee > 0L])))
  if (!is.numeric(priorities) || !length(priorities) || anyNA(priorities) ||
      any(!is.finite(priorities)) ||
      any(priorities < 1 | priorities > .Machine$integer.max |
            priorities != floor(priorities)) ||
      anyDuplicated(priorities) || priorities[1L] != 1L)
    stop("priorities must be unique positive event codes, with terminal code 1 first.",
         call. = FALSE)
  priorities <- unname(as.integer(priorities))
  if (!all(ee[ee > 0L] %in% priorities))
    stop("Every positive status code must be included in priorities.",
         call. = FALSE)
  if (is.null(labels)) {
    labels <- paste("Priority", priorities)
    labels[1L] <- "Terminal event"
  } else if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("Priority labels must be unique, nonempty names.", call. = FALSE)
  }
  d <- length(priorities)
  if (!is.numeric(thresholds) || anyNA(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0) ||
      !length(thresholds) %in% c(1L, d))
    stop("thresholds must be nonnegative, with length one or one per priority.",
         call. = FALSE)
  thresholds <- rep(as.numeric(thresholds), length.out = d)
  names(thresholds) <- as.character(priorities)
  forms <- .drw_outcome_formulas(outcome, priorities, labels)
  censoring <- .drw_formula(censoring, "censoring")
  baseline_names <- setdiff(names(data), structural)
  vars <- unique(unlist(lapply(c(forms, list(censoring)), all.vars),
                        use.names = FALSE))
  if ("." %in% vars)
    vars <- union(setdiff(vars, "."), baseline_names)
  if (length(setdiff(vars, names(data))))
    stop("Formula variables are missing from data: ",
         paste(setdiff(vars, names(data)), collapse = ", "), ".", call. = FALSE)
  if (length(intersect(vars, structural)))
    stop("Formulas must use baseline covariates, not the columns specified as ",
         "id, time, status, treatment, or cluster. Models are fitted separately ",
         "in each treatment arm.", call. = FALSE)
  for (nm in vars) {
    x <- data[[nm]]
    if (!is.atomic(x) || !is.null(dim(x)) || length(x) != nrow(data))
      stop("Formula variable ", nm, " must be an ordinary atomic column.",
           call. = FALSE)
    if (anyNA(x))
      stop("Missing values in baseline covariate ", nm,
           ". Supply a complete analysis dataset; participants are not deleted automatically.",
           call. = FALSE)
  }
  ids <- as.character(data[[id]])
  clusters <- if (is.null(cluster)) ids else as.character(data[[cluster]])
  arms_raw <- as.character(data[[treatment]])
  if (anyNA(ids) || any(!nzchar(trimws(ids))) ||
      anyNA(clusters) || any(!nzchar(trimws(clusters))) ||
      anyNA(arms_raw) || any(!nzchar(trimws(arms_raw))))
    stop("Participant IDs, cluster IDs, and treatment values must be nonmissing.",
         call. = FALSE)
  if (is.null(treatment_levels)) {
    z <- data[[treatment]]
    treatment_levels <- if (is.factor(z))
      levels(droplevels(z)) else sort(unique(arms_raw))
  }
  treatment_levels <- as.character(treatment_levels)
  if (length(treatment_levels) != 2L || anyNA(treatment_levels) ||
      anyDuplicated(treatment_levels) ||
      !setequal(treatment_levels, unique(arms_raw)))
    stop("treatment_levels must contain the two observed groups, control first ",
         "and intervention second.", call. = FALSE)
  if (any(treatment_levels == "Overall"))
    stop("The treatment label 'Overall' is reserved for censoring summaries; ",
         "rename that treatment level.", call. = FALSE)
  key <- paste(nchar(clusters), clusters, nchar(ids), ids, sep = ":")
  ord <- order(clusters, ids, tt, ee)
  keys <- unique(key[ord])
  groups <- split(seq_len(nrow(data)), factor(key, levels = keys))
  first <- vapply(groups, function(ii) ii[1L], integer(1))
  n <- length(groups)
  cluster_labels <- sort(unique(clusters[first]))
  cl <- match(clusters[first], cluster_labels)
  m <- length(cluster_labels)
  if (n < 2L || m < 2L)
    stop("At least two independent randomization units are required.",
         call. = FALSE)
  arm <- integer(n)
  stop_time <- numeric(n)
  terminal <- integer(n)
  event_times <- matrix(Inf, nrow = n, ncol = d,
                        dimnames = list(NULL, as.character(priorities)))
  repeated <- tied_terminal <- 0L
  for (i in seq_len(n)) {
    ii <- groups[[i]]
    desc <- paste0("participant ", ids[ii[1L]],
                   if (!is.null(cluster)) paste0(" in cluster ", clusters[ii[1L]]) else "")
    if (length(unique(arms_raw[ii])) != 1L)
      stop("Treatment changes within ", desc, ".", call. = FALSE)
    arm[i] <- match(arms_raw[ii[1L]], treatment_levels) - 1L
    final <- ii[ee[ii] %in% c(0L, 1L)]
    if (length(final) != 1L)
      stop("Require exactly one final record with status 0 or 1 for ", desc,
           ". A nonterminal event alone does not specify the end of follow-up.",
           call. = FALSE)
    stop_time[i] <- tt[final]
    terminal[i] <- as.integer(ee[final] == 1L)
    first[i] <- final
    if (any(tt[ii] > stop_time[i]))
      stop("Records occur after the end of follow-up for ", desc, ".",
           call. = FALSE)
    for (nm in vars) {
      if (length(unique(data[[nm]][ii])) != 1L)
        stop("Baseline covariate ", nm, " changes within ", desc, ".",
             call. = FALSE)
    }
    if (terminal[i]) event_times[i, 1L] <- stop_time[i]
    for (j in seq_len(d)[-1L]) {
      jj <- ii[ee[ii] == priorities[j]]
      if (!length(jj)) next
      if (terminal[i] && any(tt[jj] == stop_time[i])) {
        tied_terminal <- tied_terminal + sum(tt[jj] == stop_time[i])
        jj <- jj[tt[jj] < stop_time[i]]
      }
      if (!length(jj)) next
      if (length(jj) > 1L) {
        if (recurrent == "error")
          stop("More than one event with status ", priorities[j], " for ",
               desc, "; use recurrent = 'first' to analyze the first event.",
               call. = FALSE)
        repeated <- repeated + length(jj) - 1L
      }
      event_times[i, j] <- min(tt[jj])
    }
  }
  if (repeated)
    warning(repeated, " later nonterminal event record(s) were omitted; ",
            "the analysis uses the first event on each priority.",
            call. = FALSE)
  if (tied_terminal)
    warning(tied_terminal, " nonterminal event record(s) at the terminal ",
            "event time were omitted because they do not precede death.",
            call. = FALSE)
  cluster_arm <- vapply(seq_len(m), function(k) {
    a <- unique(arm[cl == k])
    if (length(a) != 1L)
      stop("Treatment is not constant within cluster ", cluster_labels[k],
           ".", call. = FALSE)
    a
  }, integer(1))
  if (!setequal(cluster_arm, 0:1))
    stop("Both treatment arms must contain an independent randomization unit.",
         call. = FALSE)
  cluster_size <- tabulate(cl, nbins = m)
  baseline <- data[first, , drop = FALSE]
  rownames(baseline) <- NULL
  model_data <- baseline[, baseline_names, drop = FALSE]
  X <- lapply(seq_len(d), function(j)
    .drw_model_matrix(forms[[j]], model_data, paste0("priority-", priorities[j], " outcome")))
  names(X) <- as.character(priorities)
  Xc <- .drw_model_matrix(censoring, model_data, "censoring")
  has_terms <- function(x) length(attr(attr(x, "terms"), "term.labels")) > 0L
  if (outcome_model %in% c("km", "nonparametric") &&
      any(vapply(X, has_terms, logical(1))))
    stop("outcome_model = '", outcome_model,
         "' requires intercept-only outcome formulas (~ 1). ",
         "Use outcome_model = 'cox' to adjust for baseline covariates.",
         call. = FALSE)
  if (censoring_model == "km" && has_terms(Xc))
    stop("censoring_model = 'km' requires censoring = ~ 1. ",
         "Use censoring_model = 'cox' to adjust for baseline covariates.",
         call. = FALSE)
  if (!is.numeric(time_probs) || !length(time_probs) || anyNA(time_probs) ||
      any(!is.finite(time_probs)) || any(time_probs <= 0 | time_probs >= 1))
    stop("time_probs must contain probabilities strictly between zero and one.",
         call. = FALSE)
  default_times <- is.null(times)
  if (default_times) {
    proposed <- as.numeric(stats::quantile(stop_time, probs = time_probs,
                                            names = FALSE, type = 7))
    times <- sort(unique(proposed[proposed > 0]))
    if (length(times) < length(proposed))
      warning("Duplicate or nonpositive default times were removed. ",
              "The fitted object reports the actual restriction times; ",
              "supply times explicitly to choose another set.", call. = FALSE)
  } else {
    if (!is.numeric(times) || !length(times) || anyNA(times) ||
        any(!is.finite(times)) || any(times <= 0))
      stop("times must contain finite, strictly positive restriction times.",
           call. = FALSE)
    times <- sort(unique(as.numeric(times)))
  }
  if (!length(times))
    stop("Positive restriction times cannot be obtained from observed follow-up.",
         call. = FALSE)
  if (!is.numeric(treatment_probability) || anyNA(treatment_probability) ||
      any(!is.finite(treatment_probability)) ||
      any(treatment_probability <= 0 | treatment_probability >= 1))
    stop("treatment_probability must be numeric and strictly between zero and one.",
         call. = FALSE)
  if (allocation == "complete") {
    if (length(treatment_probability) != 1L)
      stop("Complete randomization uses the observed fixed arm totals; ",
           "supply a scalar treatment_probability.", call. = FALSE)
    cluster_pi <- rep(mean(cluster_arm), m)
  } else if (length(treatment_probability) == 1L) {
    cluster_pi <- rep(as.numeric(treatment_probability), m)
  } else if (length(treatment_probability) == m) {
    if (!is.null(names(treatment_probability))) {
      if (anyDuplicated(names(treatment_probability)) ||
          !setequal(names(treatment_probability), cluster_labels))
        stop("Name treatment probabilities by the observed randomization units.",
             call. = FALSE)
      treatment_probability <- treatment_probability[cluster_labels]
    }
    cluster_pi <- unname(as.numeric(treatment_probability))
  } else {
    stop("treatment_probability must have length one or one value per ",
         "randomization unit.", call. = FALSE)
  }
  irt <- all(cluster_size == 1L)
  targets <- if (target == "both") {
    if (irt) "individual" else c("cluster", "individual")
  } else target
  design <- list(
    type = if (irt) "IRT" else "CRT", n = n, m = m,
    n_by_arm = setNames(tabulate(arm + 1L, nbins = 2L), treatment_levels),
    m_by_arm = setNames(tabulate(cluster_arm + 1L, nbins = 2L), treatment_levels),
    cluster_size = cluster_size, cluster_labels = cluster_labels,
    cluster_size_summary = c(min = min(cluster_size), mean = mean(cluster_size),
                             median = stats::median(cluster_size),
                             max = max(cluster_size)),
    allocation = allocation, treatment_levels = treatment_levels,
    targets = targets, times = times, default_times = default_times,
    time_probs = if (default_times) time_probs else NULL,
    time_quantile_type = if (default_times) 7L else NULL,
    id_column = id, cluster_column = cluster, time_column = time,
    status_column = status, treatment_column = treatment,
    final_status_zero = setNames(
      c(sum(terminal == 0L), sum(terminal[arm == 0L] == 0L),
        sum(terminal[arm == 1L] == 0L)),
      c("Overall", treatment_levels))
  )
  list(
    n = n, m = m, cluster = as.integer(cl), arm = as.integer(arm),
    cluster_arm = cluster_arm, cluster_size = cluster_size,
    stop = stop_time, terminal = terminal, event_times = event_times,
    baseline = baseline, X = X, Xc = Xc, treatment_levels = treatment_levels,
    times = times, priorities = priorities, labels = labels,
    thresholds = thresholds, target = target, targets = targets,
    pi = cluster_pi[cl], cluster_pi = cluster_pi,
    outcome_formulas = forms, censoring_formula = censoring,
    outcome_model = outcome_model, censoring_model = censoring_model,
    design = design,
    censoring = .drw_censoring_table(stop_time, terminal, arm, cl, times,
                                     targets, treatment_levels),
    diagnostics = list(repeated_nonterminal_records = repeated,
                       nonterminal_records_at_death = tied_terminal),
    subject_keys = keys
  )
}
