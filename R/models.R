# Nuisance models. All influence functions use the randomization unit as the
# independent unit. Positive ridge penalties are explicit user choices.
.drw_penalty <- function(control, name) {
    x <- control[[name]]
    if (is.null(x)) x <- control$ridge
    if (is.null(x)) x <- 0
    if (length(x) != 1L || !is.finite(x) || x < 0)
        stop(name, " must be a finite nonnegative number.", call. = FALSE)
    as.double(x)
}

.drw_model_design <- function(X) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    p <- ncol(X)
    if (!p) return(list(X = X, center = numeric(), scale = numeric(),
                        keep = integer(), p = 0L, names = colnames(X)))
    center <- colMeans(X)
    scale <- sqrt(colMeans(sweep(X, 2L, center)^2))
    keep <- which(is.finite(scale) & scale > 1e-12)
    if (length(keep) < p)
        warning("Constant covariates in one treatment arm have zero fitted coefficients: ",
                paste(colnames(X)[setdiff(seq_len(p), keep)], collapse = ", "),
                call. = FALSE)
    Z <- sweep(sweep(X[, keep, drop = FALSE], 2L, center[keep]),
               2L, scale[keep], "/")
    list(X = Z, center = center, scale = scale, keep = keep, p = p,
         names = colnames(X))
}

.drw_information_solve <- function(H, score, m, what) {
    q <- ncol(score)
    if (!q) return(matrix(0, m, 0L))
    H <- (H + t(H)) / 2
    if (any(!is.finite(H))) stop(what, ": nonfinite information matrix.", call. = FALSE)
    ch <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(ch))
        stop(what, ": the fitted information matrix is singular or not positive definite. ",
             "Check collinearity, simplify the model, or explicitly request a positive ridge penalty.",
             call. = FALSE)
    centered <- sweep(score, 2L, colMeans(score))
    # H is the information for the sum likelihood, while the IF is scaled by m.
    m * t(backsolve(ch, forwardsolve(t(ch), t(centered))))
}

.drw_restore_margin <- function(beta, hazard, beta_if, hazard_if, design) {
    m <- nrow(beta_if)
    K <- length(hazard)
    full_beta <- setNames(numeric(design$p), design$names)
    full_if <- matrix(0, m, design$p, dimnames = list(NULL, design$names))
    if (length(design$keep)) {
        full_beta[design$keep] <- beta / design$scale[design$keep]
        full_if[, design$keep] <- sweep(beta_if, 2L, design$scale[design$keep], "/")
    }
    factor <- exp(-sum(design$center * full_beta))
    cumulative <- cumsum(hazard) * factor
    survival <- exp(-cumulative)
    IF <- matrix(0, m, K)
    if (K) {
        cum_if <- hazard_if
        if (K > 1L) for (k in 2L:K) cum_if[, k] <- cum_if[, k] + cum_if[, k - 1L]
        hc <- factor * cum_if - outer(drop(full_if %*% design$center), cumulative)
        IF <- -sweep(hc, 2L, survival, "*")
        IF[, survival == 0] <- 0
    }
    if (any(!is.finite(full_beta)) || any(!is.finite(survival)) || any(!is.finite(IF)))
        stop("The fitted survival model is numerically nonfinite.", call. = FALSE)
    list(beta = full_beta, base_survival = survival,
         beta_if = full_if, survival_if = IF)
}

.drw_fit_cox_margin <- function(X, time, event, cluster, m, ridge,
                                control, inference, censor = FALSE, censor_before = NULL) {
    design <- .drw_model_design(X)
    Z <- design$X
    time <- as.double(time); event <- as.integer(event); cluster <- as.integer(cluster)
    p <- ncol(Z)
    penalty <- if (p) ridge / design$scale[design$keep]^2 else 0
    ordering <- if (is.null(censor_before)) {
        if (isTRUE(censor)) event == 0L else rep(FALSE, length(event))
    } else as.logical(censor_before)
    call <- function(b, scores = FALSE)
        .Call("drw_cox_likelihood", as.double(b), Z, time,
              event, cluster, as.integer(m), as.double(penalty),
              ordering, as.logical(scores), PACKAGE = "DRwinStats")
    if (!any(event)) {
        fit <- .drw_restore_margin(numeric(p), numeric(), matrix(0, m, p),
                                  matrix(0, m, 0L), design)
        fit$base_times <- numeric(); fit$converged <- TRUE; fit$method <- "cox"
        fit$diagnostics <- list(events = 0L, convergence = 0L, score = 0,
                               ridge = ridge, no_events = TRUE)
        fit$.design <- design; fit$.beta <- numeric(p); fit$.hazard <- numeric()
        return(fit)
    }
    cached_par <- NULL; cached <- NULL
    eval <- function(b) {
        if (is.null(cached_par) || !identical(b, cached_par)) {
            cached <<- call(b); cached_par <<- b
        }
        cached
    }
    if (p) {
        opt <- stats::optim(numeric(p), function(b) eval(b)$value,
                            function(b) eval(b)$gradient, method = "BFGS",
                            control = list(maxit = control$maxit, reltol = control$tol))
        beta <- opt$par
        converged <- opt$convergence == 0L
        if (!converged)
            stop("Cox model did not converge. Increase control$maxit, simplify the formula, ",
                 "or explicitly request a positive ridge penalty.", call. = FALSE)
    } else { beta <- numeric(); opt <- list(convergence = 0L); converged <- TRUE }
    ev <- call(beta, inference)
    # A residual score criterion catches unpenalized monotone likelihoods that
    # an objective-change stopping criterion can incorrectly call converged.
    if (p && (any(!is.finite(beta)) || max(abs(beta)) > 30))
        stop("The Cox coefficients indicate separation or a monotone likelihood. ",
             "Simplify the formula or explicitly request ridge regularization.", call. = FALSE)
    beta_if <- matrix(0, m, p)
    hazard_if <- matrix(0, m, length(ev$times))
    if (inference) {
        beta_if <- .drw_information_solve(ev$information, ev$score, m, "Cox model")
        hazard_if <- ev$baseline_score
        if (p && length(ev$times))
            hazard_if <- hazard_if - sweep(beta_if %*% ev$risk_mean, 2L, ev$hazard, "*")
    }
    fit <- .drw_restore_margin(beta, ev$hazard, beta_if, hazard_if, design)
    fit$base_times <- ev$times; fit$converged <- converged; fit$method <- "cox"
    fit$diagnostics <- list(events = sum(event), convergence = opt$convergence,
                           score = if (p) max(abs(ev$gradient)) else 0,
                           ridge = ridge, no_events = FALSE)
    fit$.design <- design; fit$.beta <- beta; fit$.hazard <- ev$hazard
    fit
}

.drw_fit_km_margin <- function(time, event, cluster, m, censor = FALSE,
                               inference = TRUE, censor_before = NULL) {
    time <- as.double(time); event <- as.integer(event)
    times <- sort(unique(time[event == 1L])); K <- length(times)
    survival <- numeric(K); IF <- matrix(0, m, K)
    current <- 1; current_if <- numeric(m)
    for (k in seq_along(times)) {
        t <- times[k]
        risk <- if (!is.null(censor_before)) time > t | (time == t & !censor_before) else
            if (censor) time > t | (time == t & event == 1L) else time >= t
        happened <- time == t & event == 1L
        Y <- sum(risk); D <- sum(happened)
        if (Y < D || Y == 0L) stop("Invalid Kaplan-Meier risk set.", call. = FALSE)
        q <- D / Y
        if (inference) {
            Yi <- tabulate(cluster[risk], nbins = m)
            Di <- tabulate(cluster[happened], nbins = m)
            q_if <- m * (Di - q * Yi) / Y
            current_if <- (1 - q) * current_if - current * q_if
        }
        current <- current * (1 - q)
        survival[k] <- current
        IF[, k] <- current_if
    }
    list(method = "km", beta = numeric(), base_times = times,
         base_survival = survival, beta_if = matrix(0, m, 0L),
         survival_if = IF, converged = TRUE,
         diagnostics = list(events = sum(event), no_events = !any(event),
                            convergence = 0L, score = 0, ridge = 0))
}

.drw_hessian <- function(par, gradient) {
    q <- length(par); H <- matrix(0, q, q)
    for (j in seq_len(q)) {
        h <- 1e-4 * (1 + abs(par[j]))
        plus <- minus <- par; plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
        H[, j] <- (gradient(plus) - gradient(minus)) / (2 * h)
    }
    (H + t(H)) / 2
}

.drw_fit_outcome_arm <- function(prepared, arm, method, family, control, inference) {
    ii <- which(prepared$arm == arm); n <- length(ii); m <- prepared$m
    d <- ncol(prepared$event_times)
    horizon <- if (is.null(prepared$horizon)) max(prepared$stop) else prepared$horizon
    stop_time <- pmin(prepared$stop[ii], horizon)
    observed_death <- prepared$terminal[ii] == 1L & prepared$stop[ii] <= horizon
    obs <- prepared$event_times[ii, , drop = FALSE]
    X <- lapply(prepared$X, function(x) x[ii, , drop = FALSE])
    cluster <- as.integer(prepared$cluster[ii]); ridge <- .drw_penalty(control, "outcome_ridge")
    if (method == "km" && (family != "independence" || any(vapply(X, ncol, integer(1L)) > 0L)))
        stop("Kaplan-Meier outcome margins require copula = 'independence' and outcome formulas ~ 1. ",
             "Use outcome_model = 'cox' for covariates or a joint copula fit.", call. = FALSE)
    if (method == "nonparametric" && any(vapply(X, ncol, integer(1L)) > 0L))
        stop("The nonparametric outcome model requires formulas ~ 1. Use 'cox' with covariates.", call. = FALSE)
    initialize <- vector("list", d)
    for (j in seq_len(d)) {
        event <- as.integer(is.finite(obs[, j]) & obs[, j] <= stop_time)
        time <- pmin(obs[, j], stop_time)
        censor_before <- j > 1L & observed_death & event == 0L
        initialize[[j]] <- if (method == "km")
            .drw_fit_km_margin(time, event, cluster, m, inference = inference,
                               censor_before = censor_before) else
            .drw_fit_cox_margin(X[[j]], time, event, cluster, m, ridge, control,
                               inference && family == "independence", censor_before = censor_before)
    }
    active <- sum(vapply(initialize, function(x) length(x$base_times) > 0L, logical(1L)))
    if (family == "independence" || active < 2L) {
        if (family != "independence")
            warning("The copula association is not identifiable with fewer than two observed event types; ",
                    "the fitted outcome distribution is independent.", call. = FALSE)
        if (family != "independence" && inference) {
            for (j in seq_len(d)) {
                event <- as.integer(is.finite(obs[, j]) & obs[, j] <= stop_time)
                initialize[[j]] <- .drw_fit_cox_margin(X[[j]], pmin(obs[, j], stop_time),
                    event, cluster, m, ridge, control, TRUE,
                    censor_before = j > 1L & observed_death & event == 0L)
            }
        }
        return(list(family = "independence", theta = 0,
                    beta = lapply(initialize, `[[`, "beta"),
                    base_times = lapply(initialize, `[[`, "base_times"),
                    base_survival = lapply(initialize, `[[`, "base_survival"),
                    beta_if = lapply(initialize, `[[`, "beta_if"),
                    survival_if = lapply(initialize, `[[`, "survival_if"),
                    theta_if = numeric(m), converged = TRUE, method = method,
                    diagnostics = list(family = "independence", margins = lapply(initialize, `[[`, "diagnostics"))))
    }
    designs <- lapply(initialize, `[[`, ".design")
    Z <- lapply(designs, `[[`, "X")
    times <- lapply(initialize, `[[`, "base_times")
    counts <- as.integer(lengths(times))
    p <- vapply(Z, ncol, integer(1L))
    # Endpoints with no events have unit survival and no estimable coefficients.
    for (j in which(counts == 0L)) {
        Z[[j]] <- matrix(0, n, 0L); p[j] <- 0L
        designs[[j]]$keep <- integer()
    }
    lower <- upper <- matrix(0L, n, d)
    for (j in seq_len(d)) {
        lower[, j] <- upper[, j] <- findInterval(stop_time, times[[j]])
        if (j > 1L) {
            precluded <- observed_death & !is.finite(obs[, j])
            lower[precluded, j] <- upper[precluded, j] <-
                findInterval(stop_time[precluded], times[[j]], left.open = TRUE)
        }
        ie <- which(is.finite(obs[, j]) & obs[, j] <= stop_time)
        if (length(ie)) {
            index <- match(obs[ie, j], times[[j]])
            lower[ie, j] <- index - 1L; upper[ie, j] <- index
        }
    }
    family_id <- match(family, c("independence", "gumbel", "clayton", "frank", "joe")) - 1L
    if (is.na(family_id)) stop("Unsupported copula family.", call. = FALSE)
    association_ridge <- control$association_ridge
    if (is.null(association_ridge)) association_ridge <- 0
    baseline_penalty <- control$baseline_penalty
    if (is.null(baseline_penalty)) baseline_penalty <- 0
    problem <- list(X = Z, lower = lower, upper = upper, counts = counts,
                    cluster = cluster, m = as.integer(m), family = as.integer(family_id),
                    association = TRUE, theta = 0, ridge = ridge,
                    association_ridge = as.double(association_ridge),
                    baseline_penalty = as.double(baseline_penalty),
                    ridge_weights = unlist(lapply(designs, function(z)
                        1 / z$scale[z$keep]^2), use.names = FALSE),
                    baseline_center = lapply(designs, function(z)
                        z$center[z$keep] / z$scale[z$keep]))
    # Covariates are standardized internally; the returned model uses the
    # original design-matrix scale and includes the corresponding IF change.
    beta_start <- unlist(lapply(seq_len(d), function(j)
        if (p[j]) initialize[[j]]$.beta else numeric()), use.names = FALSE)
    start_distance <- if (family == "frank") 3 else 0.5
    par <- c(beta_start, log(start_distance),
             unlist(lapply(initialize, function(x) log(pmax(x$.hazard, 1e-10))),
                    use.names = FALSE))
    association_index <- sum(p) + 1L
    q <- length(par)
    bounds_lower <- c(rep(-Inf, sum(p)), log(1e-5), rep(-35, sum(counts)))
    bounds_upper <- c(rep(Inf, sum(p)), log(100), rep(8, sum(counts)))
    eval_native <- function(x, scores = FALSE)
        .Call("drw_outcome_likelihood", as.double(x), problem,
              as.logical(scores), PACKAGE = "DRwinStats")
    cached_par <- NULL; cached <- NULL
    eval <- function(x) {
        if (is.null(cached_par) || !identical(cached_par, x)) {
            cached <<- eval_native(x); cached_par <<- x
        }
        cached
    }
    scales <- c(rep(1 / sqrt(max(1, sum(counts))), association_index), rep(1, sum(counts)))
    opt <- stats::optim(par, function(x) eval(x)$value, function(x) eval(x)$gradient,
                        method = "L-BFGS-B", lower = bounds_lower, upper = bounds_upper,
                        control = list(maxit = control$maxit, factr = 1e5,
                                       pgtol = control$tol, parscale = scales))
    result <- eval_native(opt$par, inference)
    if (opt$convergence != 0L || !isTRUE(result$valid))
        stop("The joint ", family, " outcome likelihood did not converge (", opt$message,
             "). Increase control$maxit or simplify the outcome formula; ridge is optional.",
             call. = FALSE)
    theta_shift <- if (family %in% c("gumbel", "joe")) 1 else 0
    theta <- theta_shift + exp(opt$par[association_index])
    parameter_if <- matrix(0, m, q)
    if (inference) {
        H <- .drw_hessian(opt$par, function(x) eval_native(x)$gradient)
        parameter_if <- .drw_information_solve(H, result$score, m, "Joint outcome model")
    }
    margins <- vector("list", d)
    beta_index <- 0L; hazard_index <- association_index
    for (j in seq_len(d)) {
        jb <- beta_index + seq_len(p[j]); jh <- hazard_index + seq_len(counts[j])
        hazards <- exp(opt$par[jh])
        hazard_if <- sweep(parameter_if[, jh, drop = FALSE], 2L, hazards, "*")
        margins[[j]] <- .drw_restore_margin(opt$par[jb], hazards,
                                            parameter_if[, jb, drop = FALSE], hazard_if, designs[[j]])
        beta_index <- beta_index + p[j]; hazard_index <- hazard_index + counts[j]
    }
    at_boundary <- any(abs(opt$par - bounds_lower) < 1e-5 |
                       abs(opt$par - bounds_upper) < 1e-5, na.rm = TRUE)
    if (at_boundary)
        warning("The joint outcome fit reached a numerical parameter boundary; inspect fit diagnostics.",
                call. = FALSE)
    list(family = family, theta = theta, beta = lapply(margins, `[[`, "beta"),
         base_times = times, base_survival = lapply(margins, `[[`, "base_survival"),
         beta_if = lapply(margins, `[[`, "beta_if"),
         survival_if = lapply(margins, `[[`, "survival_if"),
         theta_if = parameter_if[, association_index] * (theta - theta_shift),
         converged = TRUE, method = method,
         diagnostics = list(family = family, theta = theta, loglik = -result$value,
                            convergence = opt$convergence, score = max(abs(result$gradient * scales)),
                            ridge = ridge, association_ridge = association_ridge,
                            baseline_penalty = baseline_penalty, boundary = at_boundary,
                            events = counts, counts = opt$counts))
}

.drw_fit_models <- function(prepared, outcome_model, censoring_model, copula,
                            control, inference) {
    outcome_model <- tolower(outcome_model)
    censoring_model <- tolower(censoring_model)
    copula <- tolower(copula)
    if (copula == "independent") copula <- "independence"
    if (!outcome_model %in% c("cox", "nonparametric", "km"))
        stop("outcome_model must be 'cox', 'nonparametric', or 'km'.", call. = FALSE)
    if (!censoring_model %in% c("cox", "km"))
        stop("censoring_model must be 'cox' or 'km'.", call. = FALSE)
    if (!copula %in% c("independence", "gumbel", "clayton", "frank", "joe"))
        stop("copula must be 'independence', 'gumbel', 'clayton', 'frank', or 'joe'.", call. = FALSE)
    inference <- !(identical(inference, FALSE) || identical(inference, "none"))
    if (is.null(control$maxit)) control$maxit <- 1000L
    if (is.null(control$tol)) control$tol <- 1e-7
    if (censoring_model == "km" && ncol(prepared$Xc))
        stop("Kaplan-Meier censoring requires censoring = ~ 1; use 'cox' with covariates.", call. = FALSE)
    outcome <- censoring <- setNames(vector("list", 2L), c("arm0", "arm1"))
    for (a in 0:1) {
        ii <- which(prepared$arm == a)
        if (!length(ii)) stop("Both treatment arms must contain participants.", call. = FALSE)
        horizon <- if (is.null(prepared$horizon)) max(prepared$stop) else prepared$horizon
        fit_time <- pmin(prepared$stop[ii], horizon)
        censor_event <- as.integer(prepared$terminal[ii] == 0L & prepared$stop[ii] < horizon)
        censoring[[a + 1L]] <- if (censoring_model == "km")
            .drw_fit_km_margin(fit_time, censor_event,
                               prepared$cluster[ii], prepared$m, censor = TRUE, inference = inference) else
            .drw_fit_cox_margin(prepared$Xc[ii, , drop = FALSE], fit_time,
                                censor_event, prepared$cluster[ii], prepared$m,
                                .drw_penalty(control, "censor_ridge"), control, inference, censor = TRUE)
        outcome[[a + 1L]] <- .drw_fit_outcome_arm(prepared, a, outcome_model, copula,
                                               control, inference)
    }
    list(outcome = outcome, censoring = censoring, compute_if = inference,
         diagnostics = list(outcome = lapply(outcome, `[[`, "diagnostics"),
                            censoring = lapply(censoring, `[[`, "diagnostics")))
}
