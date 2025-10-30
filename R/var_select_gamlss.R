#' Beta regression LASSO via GAMLSS
#'
#' Uses `gamlss::ri()` (L1 penalty) in a `gamlss(dist = BE)` mean submodel to
#' select variables. The helper works on complete cases of `X`/`Y`, targets the
#' mean component, and does not yet expose offset handling.
#'
#' @inheritParams betareg_step_aic
#' @param method `"ML"` or `"GAIC"` (see `gamlss::ri`).
#' @param k Penalty multiplier for GAIC when `method = "GAIC"`.
#' @param degf Optional degrees of freedom for the L1 term.
#' @param lambda Optional penalty strength.
#' @return Named numeric vector of coefficients `(Intercept)` + `colnames(X)`,
#'   with 0 for unselected variables.
#' @seealso [gamlss::gamlss()], [gamlss::ri()], [gamlss.dist::BE()]
#' @examples
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1]); Y <- rbeta(100, Y*30, (1-Y)*30)
#' betareg_lasso_gamlss(X, Y, method = "GAIC", k = 2)
#' @export
betareg_lasso_gamlss <- function(
  X, Y, method = c("ML", "GAIC"), k = 2, degf = NULL, lambda = NULL, trace = FALSE
) {
  if (!requireNamespace("gamlss", quietly = TRUE) ||
      !requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("gamlss and gamlss.dist are required for betareg_lasso_gamlss().")
  }

  method <- match.arg(method)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  xnames <- colnames(X)
  n <- length(Y)

  # squeeze Y into (0,1)
  y <- pmin(pmax((as.numeric(Y) * (n - 1) + 0.5) / n, .Machine$double.eps), 1 - .Machine$double.eps)
  dat <- data.frame(y = y, X, check.names = TRUE)
  
  gamlss <- get("gamlss", asNamespace("gamlss"))
  ri <- get("ri", asNamespace("gamlss"))
  gamlss.control <- get("gamlss.control", asNamespace("gamlss"))
  glim.control <- get("glim.control", asNamespace("gamlss"))


  datnames <- names(dat)[-1]
  fit <- eval(substitute(gamlss(y~ri(x.vars=nnames,
                                       method=mmethod,
                                       k=kk,
                                       df=ddegf,
                                       lambda=llambda,
                                       Lp=1
  ),
  sigma.fo=~1,
  data=ddat,
  family=gamlss.dist::BE(),
  trace=ttrace
  ),
  list(nnames=datnames, mmethod=method, kk=k, ddegf=degf, llambda=lambda, ttrace=trace, ddat=dat)
  ))
  
  
  # # Use a more direct approach for the smooth term
  # formula_str <- paste("y ~", paste("ri(x.vars=c(\"", paste(xnames, collapse = "\",\""), 
  #                                   "\"), method = '", method, "', k = ", k,
  #                                   ", c.crit = 1e-06",  # CRITICAL: prevents parent environment pollution
  #                                   ", start = 10",  # CRITICAL: prevents parent environment pollution
  #                                   ", order = 0",  # CRITICAL: prevents parent environment pollution
  #                                   ", Lp = 1",     # L1 penalty for LASSO
  #                                   if (!is.null(degf)) paste(", df =", degf) else ", df = NULL",
  #                                   if (!is.null(lambda)) paste(", lambda =", lambda) else ", lambda = NULL",
  #                                   ")", sep = ""))
  

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
attr(betareg_lasso_gamlss, "fun.name") <- "betareg_lasso_gamlss"

#' Beta regression Elastic-Net via GAMLSS (gamlss.lasso)
#'
#' Uses `gamlss.lasso::gnet()` to fit ENet on the mean submodel of
#' `gamlss(dist = BE)`. The routine assumes complete cases and does not expose
#' offsets or precision-model terms.
#'
#' @inheritParams betareg_step_aic
#' @param method `"IC"` (information criterion) or `"CV"`.
#' @param ICpen Penalty for `"IC"` selection: `"BIC"`, `"AIC"`, or `"HQC"`.
#' @param alpha Elastic-net mixing (1 = LASSO, 0 = ridge).
#' @return Named numeric vector of coefficients as in [betareg_lasso_gamlss()].
#' @seealso [gamlss.lasso::gnet()], [gamlss::gamlss()], [gamlss.dist::BE()]
#' @export
betareg_enet_gamlss <- function(
  X, Y, method = c("IC", "CV"), ICpen = c("BIC", "AIC", "HQC"), alpha = 1, trace = FALSE
) {
  if (!requireNamespace("gamlss.lasso", quietly = TRUE))
    stop("Please install 'gamlss.lasso' to use betareg_enet_gamlss().")
  method <- match.arg(method); ICpen <- match.arg(ICpen)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X)
  n <- length(Y)
  
  # squeeze Y into (0,1)
  y <- pmin(pmax((as.numeric(Y) * (n - 1) + 0.5) / n, .Machine$double.eps), 1 - .Machine$double.eps)
  dat <- data.frame(y = y, X, check.names = TRUE)
  gamlss <- get("gamlss", asNamespace("gamlss"))
  gnet <- get("gnet", asNamespace("gamlss.lasso"))
  gamlss.gnet <- get("gamlss.gnet", asNamespace("gamlss.lasso"))
  gnet_control <- get("gnet.control", asNamespace("gamlss.lasso"))
  
  datnames <- names(dat)[-1]
  fit <- eval(substitute(gamlss(y~gnet(x.vars=nnames,
                                       method=mmethod,
                                       ICpen=IICpen,
                                       control=gnet_control(alpha=aalpha, 
                                                            standardize=TRUE)
                                       ),
                sigma.fo=~1, 
                data=ddat,
                family=gamlss.dist::BE(), 
                trace=ttrace
                ),
                list(nnames=datnames, mmethod=method, IICpen=ICpen, aalpha=alpha, ddat=dat, ttrace=trace)
  ))

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
attr(betareg_enet_gamlss, "fun.name") <- "betareg_enet_gamlss"
