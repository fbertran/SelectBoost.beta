#' @keywords internal
#' @aliases SelectBoost.beta-package SelectBoost.beta NULL
#'
#' @references TODO
#'
#' @author This package was written by Frédéric Bertrand.
#' Maintainer: Frédéric Bertrand <frederic.bertrand@@lecnam.net>
#'
#' @examples
#'
#' set.seed(1)
#' n <- 150; p <- 6
#' X <- matrix(rnorm(n*p), n, p); colnames(X) <- paste0("x",1:p)
#' eta <- 0.4 + X[,1] - 0.7*X[,3]
#' mu  <- plogis(eta)
#' Y   <- rbeta(n, mu*25, (1-mu)*25)
#' 
#' betareg_step_aic(X, Y)    # should return (Intercept) + x1,x3 nonzero often
#' betareg_step_bic(X, Y)
#' betareg_step_aicc(X, Y)
#'
"_PACKAGE"

# export(betareg_step_aic)
# export(betareg_step_bic)
# export(betareg_step_aicc)
# export(betareg_lasso_gamlss)
# export(betareg_enet_gamlss)
# export(betareg_glmnet)
# export(fastboost_interval)
# export(compare_selectors_single)
# export(compare_selectors_bootstrap)
# export(compare_table)
# export(plot_compare_coeff)
# export(plot_compare_freq)
# export(simulation_DATA.beta)
#' @importFrom graphics box
#' @importFrom graphics image
#' @importFrom rlang .data
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats plogis
#' @importFrom stats qbeta
#' @importFrom stats qnorm
#' @importFrom stats rbeta
#' @importFrom stats rbinom
#' @importFrom stats reshape
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats rt
#' @importFrom stats runif
#' @importFrom stats setNames
#' @importFrom utils tail
#' @useDynLib SelectBoost.beta, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
