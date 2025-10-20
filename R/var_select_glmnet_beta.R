#' Pure glmnet IRLS selector for Beta regression
#'
#' Runs an IRLS loop with Beta working responses/weights and calls
#' `glmnet` on the weighted least-squares surrogate. Supports BIC/AIC/CV
#' model choice and an optional `prestandardize` speedup.
#'
#' @inheritParams betareg_step_aic
#' @param alpha Elastic-net mixing parameter.
#' @param choose One of `"bic"`, `"aic"`, or `"cv"` to pick `lambda`.
#' @param nfolds Folds for CV when `choose = "cv"`.
#' @param n_iter Max IRLS iterations; 
#' @param tol Convergence tolerance for the IRLS parameter change (Euclidean norm
#' of the difference in `[a0, beta]`), default `1e-5`.
#' @param standardize Forwarded to `glmnet` (ignored if `prestandardize = TRUE`).
#' @param lambda Optional fixed lambda; if `NULL`, chosen by `choose`.
#' @param phi_init Initial precision (phi).
#' @param update_phi Logical; update phi inside the IRLS loop.
#' @param phi_maxit Newton steps for phi update.
#' @param prestandardize If `TRUE`, manually center/scale `X` once and
#'   disable `glmnet`'s internal standardization (speed trick).
#' @param trace Logical; print IRLS progress.
#'
#' @return Named numeric vector `(Intercept)` + `colnames(X)` with zeros for
#'   unselected variables.
#' @seealso [glmnet::glmnet()], [glmnet::cv.glmnet()]
#' @examples
#' set.seed(1); X <- matrix(rnorm(500), 100, 5); Y <- plogis(X[,1]-0.5*X[,3])
#' Y <- rbeta(100, Y*40, (1-Y)*40)
#' betareg_glmnet(X, Y, alpha = 1, choose = "bic", prestandardize = TRUE)
#' @export
betareg_glmnet <- function(
  X, Y, alpha = 1, choose = c("bic", "aic", "cv"), nfolds = 5, n_iter = 6, tol = 1e-5,
  standardize = TRUE, lambda = NULL, phi_init = 20, update_phi = TRUE, phi_maxit = 5,
  prestandardize = FALSE, trace = FALSE
) {
  choose <- match.arg(choose)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X)
  n <- nrow(X); p <- ncol(X)
  y <- .sqz01(as.numeric(Y))
  if (prestandardize) {
    xm <- colMeans(X)
    xs <- sqrt(colSums((X - rep(xm, each = n))^2) / pmax(1, n - 1))
    xs[!is.finite(xs) | xs == 0] = 1
    Xw <- (X - rep(xm, each = n)) / rep(xs, each = n)
    stan_glmnet <- FALSE
  } else {
    Xw <- X; xm <- rep(0, p); xs <- rep(1, p); stan_glmnet <- standardize
  }
  beta  <- numeric(p); a0 <- stats::qlogis(mean(y))
  eta   <- as.vector(a0 + X %*% beta)
  mu    <- stats::plogis(eta); mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  phi   <- max(phi_init, 1e-6)

  pick_lambda_ic <- function(fit, W, z, Xw) {
    pred <- drop(matrix(fit$a0, nrow = n, ncol = length(fit$lambda), byrow = TRUE) +
                 Xw %*% as.matrix(fit$beta))
    rssw <- colSums((z - pred)^2 * W)
    df   <- colSums(as.matrix(fit$beta) != 0) + 1L
    sig2 <- rssw / n
    if (choose == "bic") n * log(sig2) + log(n) * df else n * log(sig2) + 2 * df
  }

  for (iter in seq_len(n_iter)) {
    w <- beta_work_cpp(y, eta, phi); mu <- w$mu; W <- w$W; z <- w$z
    if (choose == "cv" && is.null(lambda)) {
      cv <- glmnet::cv.glmnet(x = Xw, y = z, weights = W, alpha = alpha, nfolds = nfolds,
                              family = "gaussian", standardize = stan_glmnet, intercept = TRUE)
      j <- which.min(abs(cv$glmnet.fit$lambda - cv$lambda.min))
      a0_s <- as.numeric(cv$glmnet.fit$a0[j]); b_s <- as.numeric(cv$glmnet.fit$beta[, j, drop = TRUE])
    } else {
      fit <- glmnet::glmnet(x = Xw, y = z, weights = W, alpha = alpha, family = "gaussian",
                            standardize = stan_glmnet, intercept = TRUE, lambda = lambda)
      if (is.null(lambda)) { ic <- pick_lambda_ic(fit, W, z, Xw); j <- which.min(ic) } else j <- 1L
      a0_s <- fit$a0[j]; b_s <- as.numeric(fit$beta[, j, drop = TRUE])
    }
    if (prestandardize) { a0_new <- as.numeric(a0_s - sum(b_s * (xm / xs))); b_new <- b_s / xs }
    else { a0_new <- a0_s; b_new <- b_s }
    eta_new <- as.vector(a0_new + X %*% b_new); mu_new <- stats::plogis(eta_new)
    mu_new  <- pmin(pmax(mu_new, 1e-6), 1 - 1e-6)
    if (update_phi) { phi <- beta_phi_update_cpp(phi, mu_new, y, phi_maxit); if (!is.finite(phi) || phi <= 0) phi <- phi_init }
    delta <- sqrt((a0 - a0_new)^2 + sum((beta - b_new)^2))
    if (trace) cat(sprintf("iter %d: ||\u1D6AB \u03B2||=%.3e, phi=%.3f\n", iter, delta, phi))
    a0 <- a0_new; beta <- b_new; eta <- eta_new; mu <- mu_new
    if (delta < tol) break
  }
  beta_out <- setNames(numeric(p), xnames); nz <- which(beta != 0); if (length(nz)) beta_out[nz] <- beta[nz]
  c("(Intercept)" = a0, beta_out)
}
