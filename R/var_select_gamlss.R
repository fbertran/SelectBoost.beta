#' Beta regression LASSO via GAMLSS
#'
#' Uses `gamlss::ri()` (L1 penalty) in a `gamlss(dist = BE)` mean submodel to
#' select variables.
#'
#' @inheritParams betareg_step_aic
#' @param method `"ML"` or `"GAIC"` (see `gamlss::ri`).
#' @param k Penalty multiplier for GAIC when `method = "GAIC"`.
#' @param df Optional degrees of freedom for the L1 term.
#' @param lambda Optional penalty strength.
#' @return Named numeric vector of coefficients `(Intercept)` + `colnames(X)`,
#'   with 0 for unselected variables.
#' @seealso [gamlss::gamlss()], [gamlss::ri()], [gamlss.dist::BE()]
#' @examples
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1]); Y <- rbeta(100, Y*30, (1-Y)*30)
#' betareg_lasso_gamlss(X, Y, method = "GAIC", k = 2)
#' @export
betareg_lasso_gamlss <- function(
  X, Y, method = c("ML", "GAIC"), k = 2, df = NULL, lambda = NULL, trace = FALSE
) {
  method <- match.arg(method)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X)
  n <- length(Y)
  y <- pmin(pmax((as.numeric(Y) * (n - 1) + 0.5) / n, .Machine$double.eps), 1 - .Machine$double.eps)
  dat <- data.frame(y = y, X, check.names = TRUE)
  ri_term <- gamlss::ri(x.vars = xnames, Lp = 1, df = df, lambda = lambda, method = method, k = k)
  fit <- gamlss::gamlss(formula = y ~ ri_term, sigma.fo = ~ 1,
                        family = gamlss.dist::BE(), data = dat, trace = trace)
  intercept <- unname(gamlss::coefAll(fit, what = "mu")["(Intercept)"])
  smo_list <- gamlss::getSmo(fit, "mu"); beta_hat <- tail(smo_list, 1)[[1]]$beta
  beta <- setNames(numeric(length(xnames)), xnames)
  if (!is.null(beta_hat)) beta[names(beta_hat)] <- beta_hat
  c("(Intercept)" = intercept, beta)
}

#' Beta regression Elastic-Net via GAMLSS (gamlss.lasso)
#'
#' Uses `gamlss.lasso::gnet()` to fit ENet on the mean submodel of
#' `gamlss(dist = BE)`.
#'
#' @inheritParams betareg_step_aic
#' @param method `"IC"` (information criterion) or `"CV"`.
#' @param ICpen Penalty for `"IC"` selection: `"BIC"`, `"AIC"`, or `"HQC"`.
#' @param alpha Elastic-net mixing (1 = LASSO, 0 = ridge).
#' @return Named numeric vector of coefficients as in [betareg_lasso_gamlss()].
#' @seealso [gamlss.lasso::gnet()], [gamlss::gamlss()], [gamlss.dist::BE()]
#' @examples
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1]); Y <- rbeta(100, Y*30, (1-Y)*30)
#' betareg_enet_gamlss(X, Y, method = "IC", ICpen = "BIC", alpha = 0.8)
#' @export
betareg_enet_gamlss <- function(
  X, Y, method = c("IC", "CV"), ICpen = c("BIC", "AIC", "HQC"), alpha = 1, trace = FALSE
) {
  if (!requireNamespace("gamlss.lasso", quietly = TRUE))
    stop("Please install 'gamlss.lasso' to use betareg_enet_gamlss()")
  method <- match.arg(method); ICpen <- match.arg(ICpen)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X)
  n <- length(Y)
  y <- pmin(pmax((as.numeric(Y) * (n - 1) + 0.5) / n, .Machine$double.eps), 1 - .Machine$double.eps)
  dat <- data.frame(y = y, X, check.names = TRUE)
  gnet_term <- gamlss.lasso::gnet(x.vars = xnames, method = method, ICpen = ICpen,
                                  control = gamlss.lasso::gnet.control(alpha = alpha, standardize = TRUE))
  fit <- gamlss::gamlss(formula = y ~ gnet_term, sigma.fo = ~ 1,
                        family = gamlss.dist::BE(), data = dat, trace = trace)
  intercept <- unname(gamlss::coefAll(fit, what = "mu")["(Intercept)"])
  smo_list  <- gamlss::getSmo(fit, "mu"); beta_hat <- tail(smo_list, 1)[[1]]$beta
  beta <- setNames(numeric(length(xnames)), xnames)
  if (!is.null(beta_hat)) beta[names(beta_hat)] <- beta_hat
  c("(Intercept)" = intercept, beta)
}
