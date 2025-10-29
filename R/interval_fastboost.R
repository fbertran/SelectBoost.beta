#' Interval-response stability selection (fastboost variant)
#'
#' Repeats selection on interval-valued responses by sampling a pseudo-response
#' from each interval (uniformly or midpoint), tallying variable selection
#' frequencies across `B` replicates.
#'
#' @param X Numeric matrix (n × p).
#' @param Y_low,Y_high Interval bounds in \[0,1\]. Rows with missing bounds are dropped.
#' @param func Function \code{function(X, y, ...)} returning a named coefficient
#'   vector as in the other selectors (nonselected = 0).
#' @param B Number of interval resamples.
#' @param sample `"uniform"` (default) or `"midpoint"` for drawing pseudo-responses.
#' @param version Ignored (reserved for future).
#' @param use.parallel Use `parallel::mclapply` if available.
#' @param seed Optional RNG seed.
#' @param ... Extra args forwarded to `func`.
#'
#' @return A list with:
#' \describe{
#'   \item{betas}{`B × (p+1)` matrix of coefficients over replicates.}
#'   \item{freq}{Named vector of selection frequencies for each predictor.}
#' }
#' @examples
#' # suppose you have interval data (Y_low, Y_high)
#' set.seed(1)
#' n <- 120; p <- 6
#' X <- matrix(rnorm(n*p), n, p); colnames(X) <- paste0("x",1:p)
#' mu <- plogis(X[,1] - 0.5*X[,2]); Y <- rbeta(n, mu*25, (1-mu)*25)
#' Y_low <- pmax(0, Y - 0.05); Y_high <- pmin(1, Y + 0.05)
#' fb <- fastboost_interval(X, Y_low, Y_high,
#'        func = function(X,y) betareg_glmnet(X,y, choose="bic", prestandardize=TRUE),
#'        B = 40)
#' sort(fb$freq, decreasing = TRUE)
#' @export
fastboost_interval <- function(
  X, Y_low, Y_high,
  func, B = 100, sample = c("uniform","midpoint"),
  version = "glmnet", use.parallel = FALSE, seed = NULL, ...
) {
  stopifnot(is.matrix(X))
  n <- nrow(X); p <- ncol(X)
  sample <- match.arg(sample)
  if (length(Y_low) != n || length(Y_high) != n)
    stop("Y_low and Y_high must have length nrow(X)")
  ok <- is.finite(Y_low) & is.finite(Y_high)
  X <- X[ok,,drop=FALSE]; Y_low <- Y_low[ok]; Y_high <- Y_high[ok]
  n <- nrow(X)

  if (!is.null(seed)) set.seed(seed)
  run_one <- function(b) {
    if (sample == "uniform") {
      u <- stats::runif(n); yb <- Y_low + u * (Y_high - Y_low)
    } else yb <- 0.5 * (Y_low + Y_high)
    yb <- pmin(pmax(yb, 0), 1)
    as.numeric(func(X, yb, ...))
  }
  if (use.parallel && requireNamespace("parallel", quietly = TRUE)) {
    idx <- seq_len(B)
    res <- parallel::mclapply(idx, run_one, mc.cores = max(1, parallel::detectCores() - 1L))
  } else res <- lapply(seq_len(B), run_one)
  betas <- do.call(rbind, res)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  cn <- c("(Intercept)", colnames(X)); colnames(betas) <- cn
  nz <- betas[, -1, drop = FALSE] != 0
  freq <- colMeans(nz); names(freq) <- colnames(X)
  structure(list(betas = betas, freq = freq, call = match.call()), class = "fastboost_interval")
}



#' SelectBoost workflow for interval responses
#'
#' @description
#' `sb_beta_interval()` forwards to [sb_beta()] while activating interval sampling
#' so that beta-regression SelectBoost runs can ingest lower/upper response
#' bounds directly. It mirrors [fastboost_interval()] but reuses the correlated
#' resampling pipeline of `sb_beta()`.
#'
#' @inheritParams sb_beta
#' @param sample Interval sampling scheme passed to the `interval` argument of
#'   [sb_beta()]. `"uniform"` draws a pseudo-response uniformly within each
#'   interval; `"midpoint"` always chooses the midpoint.
#' @param Y Optional point-valued response. Supply it when you wish to keep the
#'   observed mean response but still resample within `Y_low`/`Y_high` for the
#'   stability steps.
#'
#' @return See [sb_beta()]. The returned object carries the same
#'   `"sb_beta"`-class attributes describing the correlation thresholds,
#'   resampling diagnostics, selector, and number of replicates.
#'
#' @examples
#' 
#' set.seed(1)
#' sim <- simulation_DATA.beta(n = 120, p = 5, s = 2)
#' y_low <- pmax(sim$Y - 0.05, 0)
#' y_high <- pmin(sim$Y + 0.05, 1)
#' interval_fit <- sb_beta_interval(
#'   sim$X,
#'   Y_low = y_low,
#'   Y_high = y_high,
#'   B = 5,
#'   step.num = 0.4
#' )
#' attr(interval_fit, "interval")
#'
#' @export
sb_beta_interval <- function(
    X,
    Y_low,
    Y_high,
    selector = betareg_step_aic,
    sample = c("uniform", "midpoint"),
    Y = NULL,
    ...
) {
  sample <- match.arg(sample)
  sb_beta(
    X = X,
    Y = Y,
    selector = selector,
    interval = sample,
    Y_low = Y_low,
    Y_high = Y_high,
    ...
  )
}
