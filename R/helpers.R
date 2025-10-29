# Internal helpers

# Compute AIC / BIC robustly for betareg fits
.aic_betareg <- function(fit, k = 2) {
  ll <- tryCatch(stats::logLik(fit), error = function(e) NA_real_)
  if (!is.finite(ll)) return(Inf)
  kpar <- .par_count_betareg(fit)
  -2 * as.numeric(ll) + k * kpar
}

.bic_betareg <- function(fit) {
  n <- .nobs_betareg(fit)
  .aic_betareg(fit, k = log(n))
}

.sqz01 <- function(y) {
  n <- length(y)
  pmin(pmax((y * (n - 1) + 0.5) / n, .Machine$double.eps),
       1 - .Machine$double.eps)
}

.coef_betareg_full <- function(
    fit,
    mean_names,
    phi_names = NULL,
    phi_prefix = "phi|",
    phi_name_map = NULL
) {
  cf <- stats::coef(fit)
  if (is.list(cf)) {
    mean_cf <- cf[["mean"]]
    phi_cf <- cf[["precision"]]
  } else {
    mean_cf <- cf
    phi_cf <- NULL
  }
  
  mean_beta <- setNames(numeric(length(mean_names)), mean_names)
  if (length(mean_cf)) {
    common_mean <- intersect(names(mean_cf), mean_names)
    if (length(common_mean)) {
      mean_beta[common_mean] <- mean_cf[common_mean]
    }
  }
  mean_intercept <- if ("(Intercept)" %in% names(mean_cf)) unname(mean_cf["(Intercept)"]) else 0
  
  phi_beta <- NULL
  phi_labels <- NULL
  if (!is.null(phi_names) && length(phi_names)) {
    phi_labels <- paste0(phi_prefix, phi_names)
    if (!is.null(phi_name_map)) {
      phi_labels <- paste0(phi_prefix, unname(phi_name_map[phi_names]))
    }
    phi_beta <- setNames(numeric(length(phi_names)), phi_labels)
    if (length(phi_cf)) {
      # betareg names the precision component either "precision" or "phi"
      phi_coef_vals <- phi_cf
      common_phi <- intersect(names(phi_coef_vals), phi_names)
      if (length(common_phi)) {
        idx <- match(common_phi, phi_names)
        phi_beta[idx] <- phi_coef_vals[common_phi]
      }
    }
  }
  
  phi_intercept <- NULL
  if (!is.null(phi_cf) && length(phi_cf)) {
    name <- paste0(phi_prefix, "(Intercept)")
    val <- if ("(Intercept)" %in% names(phi_cf)) unname(phi_cf["(Intercept)"]) else 0
    phi_intercept <- setNames(val, name)
  }
  
  c(`(Intercept)` = mean_intercept, mean_beta, phi_intercept, phi_beta)
}

.scope_beta <- function(xnames) {
  full <- as.formula(paste("y ~", paste(xnames, collapse = " + ")))
  list(lower = y ~ 1, upper = full)
}

.par_count_betareg <- function(fit) {
  cf <- coef(fit)
  if (is.list(cf)) sum(lengths(cf)) else length(cf)
}

.nobs_betareg <- function(fit) {
  n <- tryCatch(stats::nobs(fit), error = function(e) NA_integer_)
  if (is.na(n)) n <- length(stats::fitted(fit))
  n
}

.aicc_betareg <- function(fit) {
  k <- .par_count_betareg(fit)
  n <- .nobs_betareg(fit)
  aic <- stats::AIC(fit)
  denom <- n - k - 1
  if (denom <= 0) return(Inf)
  aic + (2 * k * (k + 1)) / denom
}

.beta_formula_from_vars <- function(mean_vars, phi_vars = NULL) {
  mean_part <- if (length(mean_vars) == 0) "1" else paste(mean_vars, collapse = " + ")
  if (is.null(phi_vars)) {
    as.formula(paste0("y ~ ", mean_part))
  } else {
    phi_part <- if (length(phi_vars) == 0) "1" else paste(phi_vars, collapse = " + ")
    as.formula(paste0("y ~ ", mean_part, " | ", phi_part))
  }
}

.shorten_colnames <- function(X, maxlen = 30) {
  xn <- colnames(X)
  if (is.null(xn)) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    return(list(X = X, map = setNames(colnames(X), colnames(X))))
  }
  xn_short <- abbreviate(xn, minlength = min(maxlen, max(nchar(xn))))
  xn_short <- make.names(xn_short, unique = TRUE)
  map <- setNames(xn, xn_short)  # short -> original
  colnames(X) <- xn_short
  list(X = X, map = map)
}


.sb_future_warning_emitted <- FALSE

.sb_parallel_map <- function(indices, FUN, use.parallel = FALSE, seed = NULL) {
  if (!length(indices)) {
    return(list())
  }
  if (!isTRUE(use.parallel) || length(indices) < 2L) {
    return(lapply(indices, FUN))
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    if (!.sb_future_warning_emitted) {
      warning(
        "`use.parallel = TRUE` requires the `future.apply` package; falling back to sequential computation.",
        call. = FALSE
      )
      .sb_future_warning_emitted <<- TRUE
    }
    return(lapply(indices, FUN))
  }
  seed_spec <- if (is.null(seed)) TRUE else seed
  future.apply::future_lapply(indices, FUN, future.seed = seed_spec)
}
