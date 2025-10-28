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

.coef_glmnet_style <- function(fit, xnames) {
  cf <- stats::coef(fit)
  mean_cf <- if (is.list(cf)) cf$mean else cf
  beta <- setNames(numeric(length(xnames)), xnames)
  common <- intersect(names(mean_cf), xnames)
  if (length(common)) beta[common] <- mean_cf[common]
  intercept <- if ("(Intercept)" %in% names(mean_cf)) unname(mean_cf["(Intercept)"]) else 0
  c("(Intercept)" = intercept, beta)
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

.beta_formula_from_vars <- function(vars) {
  if (length(vars) == 0) as.formula("y ~ 1")
  else as.formula(paste0("y ~ ", paste(vars, collapse = " + ")))
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
