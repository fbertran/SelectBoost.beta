#' Stepwise Beta regression by AIC
#'
#' Fits a mean-submodel Beta regression with `betareg` and performs
#' stepwise subset selection using AIC.
#'
#' @param X Numeric matrix (n × p) of predictors.
#' @param Y Numeric response in (0,1). Values are squeezed to (0,1) internally.
#' @param direction Stepwise direction: `"both"`, `"forward"`, or `"backward"`.
#' @param link Link for the mean submodel (passed to `betareg`), default `"logit"`.
#' @param link.phi Link for precision parameter, default `"log"`.
#' @param type Likelihood type for `betareg`, e.g. `"ML"`.
#' @param trace Logical; print stepwise trace.
#' @param max_steps Integer; maximum number of greedy steps (default `p`).
#' @param epsilon Numeric; minimum AIC improvement required to accept a move (default `1e-8`).
#'
#' @return Named numeric vector of length `p+1` with `(Intercept)` and one
#'   coefficient per column of `X`. Nonselected variables have 0.
#' @seealso [betareg::betareg()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200), 100, 2); Y <- plogis(0.5 + X[,1]-X[,2]); Y <- rbeta(100, Y*20, (1-Y)*20)
#' betareg_step_aic(X, Y)
#' @export
betareg_step_aic <- function(
    X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
    trace = FALSE, max_steps = NULL, epsilon = 1e-8
) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X); if (is.null(max_steps)) max_steps <- ncol(X)
  
  y <- .sqz01(as.numeric(Y))
  dat <- data.frame(y = y, X, check.names = TRUE)
  
  current_vars <- character(0)
  current_fit  <- betareg::betareg(.beta_formula_from_vars(current_vars), data = dat,
                                   link = link, link.phi = link.phi, type = type)
  current_aic  <- .aic_betareg(current_fit, k = 2)
  
  can_forward  <- direction %in% c("both","forward")
  can_backward <- direction %in% c("both","backward")
  
  step <- 0L
  repeat {
    step <- step + 1L; if (step > max_steps) break
    candidates <- list(); cand_info <- list()
    
    if (can_forward) {
      to_add <- setdiff(xnames, current_vars)
      for (v in to_add) {
        vars <- sort(c(current_vars, v))
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (!inherits(fit, "try-error")) { candidates[[paste0("+", v)]] <- fit; cand_info[[paste0("+", v)]] <- vars }
      }
    }
    if (can_backward && length(current_vars)) {
      for (v in current_vars) {
        vars <- setdiff(current_vars, v)
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (!inherits(fit, "try-error")) { candidates[[paste0("-", v)]] <- fit; cand_info[[paste0("-", v)]] <- vars }
      }
    }
    if (!length(candidates)) break
    
    aics <- vapply(candidates, .aic_betareg, numeric(1), k = 2)
    j <- which.min(aics); key <- names(aics)[j]
    best_fit <- candidates[[key]]; best_vars <- cand_info[[key]]; best_aic <- aics[j]
    
    if (trace) message(sprintf("step %d: best %s, AIC: %.4f -> %.4f",
                               step, key, current_aic, best_aic))
    if (is.finite(best_aic) && (current_aic - best_aic) > epsilon) {
      current_fit <- best_fit; current_vars <- best_vars; current_aic <- best_aic
    } else break
  }
  
  .coef_glmnet_style(current_fit, xnames)
}

#' Stepwise Beta regression by BIC
#' @inheritParams betareg_step_aic
#' @return See [betareg_step_aic()].
#' @examples
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1]); Y <- rbeta(100, Y*30, (1-Y)*30)
#' betareg_step_bic(X, Y)
#' @export
betareg_step_bic <- function(
    X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
    trace = FALSE, max_steps = NULL, epsilon = 1e-8
) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X); if (is.null(max_steps)) max_steps <- ncol(X)
  
  y <- .sqz01(as.numeric(Y))
  dat <- data.frame(y = y, X, check.names = TRUE)
  
  current_vars <- character(0)
  current_fit  <- betareg::betareg(.beta_formula_from_vars(current_vars), data = dat,
                                   link = link, link.phi = link.phi, type = type)
  current_bic  <- .bic_betareg(current_fit)
  
  can_forward  <- direction %in% c("both","forward")
  can_backward <- direction %in% c("both","backward")
  
  step <- 0L
  repeat {
    step <- step + 1L; if (step > max_steps) break
    candidates <- list(); cand_info <- list()
    
    if (can_forward) {
      to_add <- setdiff(xnames, current_vars)
      for (v in to_add) {
        vars <- sort(c(current_vars, v))
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (!inherits(fit, "try-error")) { candidates[[paste0("+", v)]] <- fit; cand_info[[paste0("+", v)]] <- vars }
      }
    }
    if (can_backward && length(current_vars)) {
      for (v in current_vars) {
        vars <- setdiff(current_vars, v)
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (!inherits(fit, "try-error")) { candidates[[paste0("-", v)]] <- fit; cand_info[[paste0("-", v)]] <- vars }
      }
    }
    if (!length(candidates)) break
    
    bics <- vapply(candidates, .bic_betareg, numeric(1))
    j <- which.min(bics); key <- names(bics)[j]
    best_fit <- candidates[[key]]; best_vars <- cand_info[[key]]; best_bic <- bics[j]
    
    if (trace) message(sprintf("step %d: best %s, BIC: %.4f -> %.4f",
                               step, key, current_bic, best_bic))
    if (is.finite(best_bic) && (current_bic - best_bic) > epsilon) {
      current_fit <- best_fit; current_vars <- best_vars; current_bic <- best_bic
    } else break
  }
  
  .coef_glmnet_style(current_fit, xnames)
}


#' Stepwise Beta regression by AICc (finite-sample corrected AIC)
#'
#' Greedy forward/backward search minimizing AICc computed on `betareg` fits.
#'
#' @inheritParams betareg_step_aic
#' @param max_steps Maximum number of greedy steps (default `p`).
#' @param epsilon Minimal AICc improvement to accept a move.
#' @return See [betareg_step_aic()].
#' @examples
#' set.seed(1); X <- matrix(rnorm(400), 100, 4); Y <- plogis(X[,1]+0.5*X[,2])
#' Y <- rbeta(100, Y*25, (1-Y)*25); betareg_step_aicc(X, Y)
#' @export
betareg_step_aicc <- function(
  X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
  trace = FALSE, max_steps = NULL, epsilon = 1e-8
) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  xnames <- colnames(X); if (is.null(max_steps)) max_steps <- ncol(X)
  y <- .sqz01(as.numeric(Y))
  dat <- data.frame(y = y, X, check.names = TRUE)
  current_vars <- character(0)
  current_fit <- betareg::betareg(.beta_formula_from_vars(current_vars), data = dat,
                                  link = link, link.phi = link.phi, type = type)
  current_aicc <- .aicc_betareg(current_fit)
  can_forward  <- direction %in% c("both", "forward")
  can_backward <- direction %in% c("both", "backward")
  step <- 0L
  repeat {
    step <- step + 1L; if (step > max_steps) break
    candidates <- list(); cand_info <- list()
    if (can_forward) {
      to_add <- setdiff(xnames, current_vars)
      for (v in to_add) {
        vars <- sort(c(current_vars, v))
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (inherits(fit, "try-error")) next
        candidates[[paste0("+", v)]] <- fit; cand_info[[paste0("+", v)]] <- vars
      }
    }
    if (can_backward && length(current_vars)) {
      for (v in current_vars) {
        vars <- setdiff(current_vars, v)
        fml  <- .beta_formula_from_vars(vars)
        fit  <- try(betareg::betareg(fml, data = dat, link = link, link.phi = link.phi, type = type), silent = TRUE)
        if (inherits(fit, "try-error")) next
        candidates[[paste0("-", v)]] <- fit; cand_info[[paste0("-", v)]] <- vars
      }
    }
    if (length(candidates) == 0) break
    aiccs <- vapply(candidates, .aicc_betareg, numeric(1))
    best_idx <- which.min(aiccs); best_key <- names(aiccs)[best_idx]
    best_fit <- candidates[[best_key]]; best_vars <- cand_info[[best_key]]; best_aicc <- aiccs[best_idx]
    if (is.finite(best_aicc) && (current_aicc - best_aicc) > epsilon) {
      current_fit  <- best_fit; current_vars <- best_vars; current_aicc <- best_aicc
    } else break
  }
  .coef_glmnet_style(current_fit, xnames)
}
