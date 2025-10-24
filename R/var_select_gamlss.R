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
  gamlss <- get("gamlss", asNamespace("gamlss"))
  ri <- get("ri", asNamespace("gamlss"))
  smooth_term <- substitute(
    ri(
      x.vars = XN,
      Lp = 1, df = DF, lambda = LAMBDA,
      method = METHOD, k = K
    ),
    list(
      XN = xnames, DF = df, LAMBDA = lambda,
      METHOD = method, K = k
    )
  )
  call_env <- list2env(list(dat = dat), parent = environment())
  fit <- eval(
    quote(gamlss(
      formula = as.formula(
        substitute(y ~ SMOOTH, list(SMOOTH = smooth_term)),
        env = environment()
      ),
      sigma.fo = ~ 1,
      family = gamlss.dist::BE(), data = dat, trace = trace
    )),
    envir = call_env
  )
  mu_coef <- gamlss::coefAll(fit, what = "mu")[[1]]
  intercept <- unname(mu_coef["(Intercept)"])
  smo_list <- gamlss::getSmo(fit, "mu")
  if (length(smo_list) == 0) {
    beta_hat <- NULL
  } else if (!is.null(smo_list$coef)) {
    beta_hat <- drop(smo_list$coef)
  } else {
    last_smo <- tail(smo_list, 1)[[1]]
    beta_hat <- if (!is.null(last_smo$beta)) {
      last_smo$beta
    } else if (!is.null(last_smo$coef)) {
      drop(last_smo$coef)
    } else {
      NULL
    }
  }
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
  gamlss <- get("gamlss", asNamespace("gamlss"))
  gnet <- get("gnet", asNamespace("gamlss.lasso"))
  gnet_control <- get("gnet.control", asNamespace("gamlss.lasso"))
  smooth_term <- substitute(
    gnet(
      x.vars = XN, method = METHOD, ICpen = ICPEN,
      control = gnet_control(alpha = ALPHA, standardize = TRUE)    ),
    list(
      XN = xnames, METHOD = method, ICPEN = ICpen,
      ALPHA = alpha
    )
  )
  call_env <- list2env(list(dat = dat), parent = environment())
  fit <- eval(
    quote(gamlss(
      formula = as.formula(
        substitute(y ~ SMOOTH, list(SMOOTH = smooth_term)),
        env = environment()
      ),
      sigma.fo = ~ 1,
      family = gamlss.dist::BE(), data = dat, trace = trace
    )),
    envir = call_env
  )
  mu_coef <- gamlss::coefAll(fit, what = "mu")[[1]]
  intercept <- unname(mu_coef["(Intercept)"])
  smo_list <- gamlss::getSmo(fit, "mu")
  if (length(smo_list) == 0) {
    beta_hat <- NULL
  } else if (!is.null(smo_list$coef)) {
    beta_hat <- drop(smo_list$coef)
  } else {
    last_smo <- tail(smo_list, 1)[[1]]
    beta_hat <- if (!is.null(last_smo$beta)) {
      last_smo$beta
    } else if (!is.null(last_smo$coef)) {
      drop(last_smo$coef)
    } else {
      NULL
    }
  }
  beta <- setNames(numeric(length(xnames)), xnames)
  if (!is.null(beta_hat)) beta[names(beta_hat)] <- beta_hat
  c("(Intercept)" = intercept, beta)
}
