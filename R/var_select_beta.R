#' Stepwise Beta regression by AIC
#'
#' Fits a Beta regression with optional joint selection of the mean and
#' precision (phi) submodels using `betareg::betareg()`. The routine performs
#' greedy forward/backward search using the requested information criterion and
#' returns coefficients aligned with the supplied design matrix. The selectors
#' currently target the mean submodel only, require complete cases, and do not
#' expose offsets. Observation weights are passed through to `betareg()` when
#' provided.
#'
#' @param X Numeric matrix (n × p) of mean-submodel predictors.
#' @param Y Numeric response in (0,1). Values are squeezed to (0,1) internally.
#' @param direction Stepwise direction for the mean submodel: `"both"`,
#'   `"forward"`, or `"backward"`.
#' @param link Link for the mean submodel (passed to `betareg`). Default
#'   `"logit"`.
#' @param link.phi Link for precision parameter. Default `"log"`.
#' @param type Likelihood type for `betareg`, e.g. `"ML"`.
#' @param trace Logical; print stepwise trace.
#' @param max_steps Integer; maximum number of greedy steps (default `p`).
#' @param epsilon Numeric; minimum improvement required to accept a move
#'   (default `1e-8`).
#' @param X_phi Optional matrix of candidate predictors for the precision (phi)
#'   submodel. When `direction_phi` enables precision updates and `X_phi` is
#'   `NULL`, the function reuses `X`.
#' @param direction_phi Stepwise direction for the precision submodel.
#'   Defaults to `"none"` (no phi selection). Supported values mirror
#'   `direction`.
#' @param weights Optional non-negative observation weights passed to
#'   `betareg()`.
#'
#' @return Named numeric vector of length `p_mean + p_phi + 1` containing the
#'   intercept, mean coefficients, phi-intercept (prefixed by `"phi|"`), and
#'   phi coefficients (also prefixed by `"phi|"`). Non-selected variables have
#'   coefficient 0.
#' @seealso [betareg::betareg()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200), 100, 2);
#' Y <- plogis(0.5 + X[,1]-X[,2]);
#' betareg_step_aic(X, Y)
#' Y <- rbeta(100, Y*20, (1-Y)*20)
#' betareg_step_aic(X, Y)
#' @export
betareg_step_aic <- function(
    X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
    trace = FALSE, max_steps = NULL, epsilon = 1e-8,
    X_phi = NULL, direction_phi = c("none", "both", "forward", "backward"),
    weights = NULL
) {
  direction_phi <- match.arg(direction_phi)
  .betareg_step_engine(
    X = X,
    Y = Y,
    criterion_fun = function(fit) .aic_betareg(fit, k = 2),
    criterion_name = "AIC",
    direction = direction,
    link = link,
    link.phi = link.phi,
    type = type,
    trace = trace,
    max_steps = max_steps,
    epsilon = epsilon,
    X_phi = X_phi,
    direction_phi = direction_phi,
    weights = weights
  )
}
attr(betareg_step_aic, "fun.name") <- "betareg_step_aic"

#' Stepwise Beta regression by BIC
#' @inheritParams betareg_step_aic
#' @return See [betareg_step_aic()].
#' @examples
#' set.seed(1); X <- matrix(rnorm(300), 100, 3);
#' Y <- plogis(X[,1]);
#' betareg_step_bic(X, Y)
#' Y <- rbeta(100, Y*30, (1-Y)*30)
#' betareg_step_bic(X, Y)
#' @export
betareg_step_bic <- function(
    X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
    trace = FALSE, max_steps = NULL, epsilon = 1e-8,
    X_phi = NULL, direction_phi = c("none", "both", "forward", "backward"),
    weights = NULL
) {
  direction_phi <- match.arg(direction_phi)
  .betareg_step_engine(
    X = X,
    Y = Y,
    criterion_fun = .bic_betareg,
    criterion_name = "BIC",
    direction = direction,
    link = link,
    link.phi = link.phi,
    type = type,
    trace = trace,
    max_steps = max_steps,
    epsilon = epsilon,
    X_phi = X_phi,
    direction_phi = direction_phi,
    weights = weights
  )
}
attr(betareg_step_bic, "fun.name") <- "betareg_step_bic"


#' Stepwise Beta regression by AICc (finite-sample corrected AIC)
#'
#' Greedy forward/backward search minimizing AICc computed on `betareg` fits with
#' optional precision-submodel selection and observation weights.
#'
#' @inheritParams betareg_step_aic
#' @param max_steps Maximum number of greedy steps (default `p`).
#' @param epsilon Minimal AICc improvement to accept a move.
#' @return See [betareg_step_aic()].
#' @examples
#' set.seed(1);
#' X <- matrix(rnorm(400), 100, 4);
#' Y <- plogis(X[,1]+0.5*X[,2])
#' betareg_step_aicc(X, Y)
#' Y <- rbeta(100, Y*25, (1-Y)*25);
#' betareg_step_aicc(X, Y)
#' @export
betareg_step_aicc <- function(
  X, Y, direction = "both", link = "logit", link.phi = "log", type = "ML",
  trace = FALSE, max_steps = NULL, epsilon = 1e-8,
  X_phi = NULL, direction_phi = c("none", "both", "forward", "backward"),
  weights = NULL
  ) {
    direction_phi <- match.arg(direction_phi)
    .betareg_step_engine(
      X = X,
      Y = Y,
      criterion_fun = .aicc_betareg,
      criterion_name = "AICc",
      direction = direction,
      link = link,
      link.phi = link.phi,
      type = type,
      trace = trace,
      max_steps = max_steps,
      epsilon = epsilon,
      X_phi = X_phi,
      direction_phi = direction_phi,
      weights = weights
    )
  }
attr(betareg_step_aicc, "fun.name") <- "betareg_step_aicc"


# Internal shared engine -------------------------------------------------------

.betareg_step_engine <- function(
    X,
    Y,
    criterion_fun,
    criterion_name,
    direction,
    link,
    link.phi,
    type,
    trace,
    max_steps,
    epsilon,
    X_phi,
    direction_phi,
    weights
) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) storage.mode(X) <- "double"
  if (!is.matrix(X)) stop("`X` must be coercible to a numeric matrix.")
  n <- nrow(X)
  if (n == 0L) stop("`X` must have at least one row.")
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  mean_names <- colnames(X)
  
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != n) {
      stop("`weights` must have length equal to nrow(X).")
    }
    if (any(!is.finite(weights) | weights < 0)) {
      stop("`weights` must be non-negative finite numbers.")
    }
  }
#  obs_weights <- weights
  
  if (direction_phi == "none") {
    X_phi <- NULL
  } else if (is.null(X_phi)) {
    X_phi <- X
  }
  
  phi_names <- character(0)
  phi_data_names <- character(0)
  phi_map <- NULL
  if (!is.null(X_phi)) {
    if (!is.matrix(X_phi)) X_phi <- as.matrix(X_phi)
    if (nrow(X_phi) != n) {
      stop("`X_phi` must have the same number of rows as `X`.")
    }
    if (is.null(colnames(X_phi))) {
      colnames(X_phi) <- paste0("Z", seq_len(ncol(X_phi)))
    }
    colnames(X_phi) <- make.names(colnames(X_phi), unique = TRUE)
    phi_names <- colnames(X_phi)
    phi_data_names <- paste0(".phi_", seq_along(phi_names))
    phi_map <- setNames(phi_names, phi_data_names)
  }
  
  y <- .sqz01(as.numeric(Y))
  dat <- data.frame(y = y, X, check.names = TRUE)
  
  if (length(phi_data_names)) {
    phi_df <- as.data.frame(X_phi, stringsAsFactors = FALSE)
    names(phi_df) <- phi_data_names
    dat <- cbind(dat, phi_df)
  }
  
  mean_candidates <- mean_names
  phi_candidates <- phi_data_names
  
  mean_vars <- character(0)
  phi_vars <- if (length(phi_candidates)) character(0) else NULL
  
  build_formula <- function(mean_set, phi_set) {
    if (!length(phi_candidates)) phi_set <- NULL
    .beta_formula_from_vars(mean_set, phi_set)
  }
  
  betareg_base_args <- list(
    data = dat,
    link = link,
    link.phi = link.phi,
    type = type
  )
  if (!is.null(weights)) {
    betareg_base_args$weights <- weights
  }
  
  fit_betareg <- function(mean_set, phi_set) {
    args_br <- c(list(formula = build_formula(mean_set, phi_set)), betareg_base_args)
    do.call(betareg::betareg, args_br)
  }
  
  current_fit <- try(
    fit_betareg(mean_vars, phi_vars),
    silent = TRUE
  )
  if (inherits(current_fit, "try-error")) {
    stop("Initial betareg fit failed: ", current_fit)
  }
  current_score <- criterion_fun(current_fit)
  
  can_forward <- direction %in% c("both", "forward")
  can_backward <- direction %in% c("both", "backward")
  can_forward_phi <- direction_phi %in% c("both", "forward") && length(phi_candidates)
  can_backward_phi <- direction_phi %in% c("both", "backward") && length(phi_candidates)
  
  if (is.null(max_steps)) {
    max_steps <- length(mean_candidates) + length(phi_candidates)
  }
  max_steps <- as.integer(max_steps)
  if (max_steps < 0L) max_steps <- 0L
  
  step <- 0L
  repeat {
    step <- step + 1L
    if (step > max_steps) break
    
    candidates <- list()
    candidate_info <- list()
    
    if (can_forward) {      
      to_add <- setdiff(mean_candidates, mean_vars)
      for (v in to_add) {
        m_vars <- sort(c(mean_vars, v))
        fit <- try(
          fit_betareg(m_vars, phi_vars),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) {
          key <- paste0("+mean:", v)
          candidates[[key]] <- fit
          candidate_info[[key]] <- list(mean = m_vars, phi = phi_vars)
        }
      }
    }
    
    if (can_backward && length(mean_vars)) {
      for (v in mean_vars) {
        m_vars <- setdiff(mean_vars, v)
        fit <- try(
          fit_betareg(m_vars, phi_vars),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) {
          key <- paste0("-mean:", v)
          candidates[[key]] <- fit
          candidate_info[[key]] <- list(mean = m_vars, phi = phi_vars)
        }
      }
    }
    
    if (can_forward_phi) {
      to_add_phi <- setdiff(phi_candidates, phi_vars)
      for (v in to_add_phi) {
        p_vars <- sort(c(phi_vars, v))
        fit <- try(
          fit_betareg(mean_vars, p_vars),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) {
          key <- paste0("+phi:", v)
          candidates[[key]] <- fit
          candidate_info[[key]] <- list(mean = mean_vars, phi = p_vars)
        }
      }
    }
    
    if (can_backward_phi && length(phi_vars)) {
      for (v in phi_vars) {
        p_vars <- setdiff(phi_vars, v)
        fit <- try(
          fit_betareg(mean_vars, p_vars),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) {
          key <- paste0("-phi:", v)
          candidates[[key]] <- fit
          candidate_info[[key]] <- list(mean = mean_vars, phi = p_vars)
        }
      }
    }
    
    
    if (!length(candidates)) break
    
    scores <- vapply(candidates, criterion_fun, numeric(1))
    best_idx <- which.min(scores)
    key <- names(scores)[best_idx]
    best_fit <- candidates[[key]]
    best_info <- candidate_info[[key]]
    best_score <- scores[best_idx]
    
    if (trace) {
      message(
        sprintf(
          "step %d: best %s, %s: %.4f -> %.4f",
          step,
          key,
          criterion_name,
          current_score,
          best_score
        )
      )
    }
    
    if (is.finite(best_score) && (current_score - best_score) > epsilon) {
      current_fit <- best_fit
      mean_vars <- best_info$mean
      phi_vars <- best_info$phi
      current_score <- best_score
    } else {
      break
    }
  }
  
  .coef_betareg_full(
    current_fit,
    mean_names = mean_names,
    phi_names = phi_vars,
    phi_name_map = phi_map
  )
}

  
  
  
  
  