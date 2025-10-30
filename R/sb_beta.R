#' Core helpers for SelectBoost-style beta regression
#'
#' These helpers expose the individual stages of the SelectBoost workflow so
#' that beta-regression selectors can be combined with correlation-aware
#' resampling directly from `SelectBoost.beta`. They normalise the design matrix,
#' derive correlation structures, form groups of correlated predictors, generate
#' Gaussian surrogates that mimic the observed dependency structure, and apply a
#' user-provided selector on each resampled design.
#'
#' @param X Numeric matrix of predictors.
#' @param center Optional centering vector recycled to the number of columns.
#'   Defaults to the column means of `X`.
#' @param scale Optional scaling vector recycled to the number of columns.
#'   Defaults to the column-wise \eqn{\ell_2} norms of the centred matrix.
#' @param eps Small positive constant used when normalising columns.
#' @return `sb_normalize()` returns a centred, \eqn{\ell_2}-scaled copy of `X`.
#' @examples
#' sb_normalize(matrix(rnorm(20), 5))
#' @export
sb_normalize <- function(X, center = NULL, scale = NULL, eps = 1e-8) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(X)) {
    storage.mode(X) <- "double"
  }
  p <- ncol(X)
  if (p == 0L) {
    attr(X, "center") <- numeric(0)
    attr(X, "scale") <- numeric(0)
    return(X)
  }
  if (is.null(center)) {
    center_vec <- colMeans(X, na.rm = TRUE)
  } else {
    center_vec <- rep_len(as.numeric(center), p)
  }
  if (!is.null(colnames(X))) {
    names(center_vec) <- colnames(X)
  }
  X_centered <- sweep(X, 2L, center_vec, FUN = "-", check.margin = FALSE)
  col_l2 <- function(v) {
    v <- v[is.finite(v)]
    norm <- sqrt(sum(v * v))
    if (!is.finite(norm) || norm < eps) 1 else norm
  }
  if (is.null(scale)) {
    scale_vec <- vapply(seq_len(p), function(j) col_l2(X_centered[, j]), numeric(1))
  } else {
    scale_vec <- rep_len(as.numeric(scale), p)
    scale_vec[!is.finite(scale_vec) | scale_vec < eps] <- 1
  }
  if (!is.null(colnames(X))) {
    names(scale_vec) <- colnames(X)
  }
  X_scaled <- sweep(X_centered, 2L, scale_vec, FUN = "/", check.margin = FALSE)
  attr(X_scaled, "center") <- center_vec
  attr(X_scaled, "scale") <- scale_vec
  X_scaled
}

#' @rdname sb_normalize
#' @param corrfunc Function or character string used to compute pairwise
#'   associations. Defaults to `"cor"`.
#' @return `sb_compute_corr()` returns the association matrix.
#' @export
sb_compute_corr <- function(X, corrfunc = "cor") {
  if (is.character(corrfunc)) {
    fun <- get(corrfunc, mode = "function", envir = parent.frame())
    fun(X)
  } else if (is.function(corrfunc)) {
    corrfunc(X)
  } else {
    stop("`corrfunc` must be a function or the name of one.")
  }
}

#' @rdname sb_normalize
#' @param corr_mat Numeric matrix of associations.
#' @param c0 Threshold applied to the absolute correlations.
#' @return `sb_group_variables()` returns a list of integer vectors, one per
#'   variable, describing the correlated group it belongs to.
#' @export
sb_group_variables <- function(corr_mat, c0) {
  cm <- abs(corr_mat)
  cm[cm < c0] <- 0
  diag(cm) <- 1
  cm[cm != 0] <- 1
  dete <- function(x) which(x == 1)
  res <- apply(cm, 2L, dete)
  if (is.vector(res)) {
    res <- lapply(res, function(x) x)
  }
  if (is.matrix(res)) {
    cn <- colnames(res)
    res <- lapply(seq_len(ncol(res)), function(i) res[, i])
    names(res) <- cn
  }
  attr(res, "type") <- "normal"
  res
}

# Internal: add jitter and ensure covariance is positive definite
.sb_cov_adjust <- function(Sigma, jitter = 1e-6) {
  Sigma[!is.finite(Sigma)] <- 0
  diag(Sigma) <- diag(Sigma) + jitter
  if (!.is_positive_definite(Sigma)) {
    eig <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
    eig$values[eig$values < jitter] <- jitter
    Sigma <- eig$vectors %*% diag(eig$values, nrow = nrow(Sigma)) %*% t(eig$vectors)
  }
  Sigma
}

.is_positive_definite <- function(Sigma) {
  if (!is.matrix(Sigma)) return(FALSE)
  if (nrow(Sigma) != ncol(Sigma)) return(FALSE)
  res <- try(chol(Sigma), silent = TRUE)
  !inherits(res, "try-error")
}

# Internal: sample correlated columns for a single group
.sb_sample_group <- function(X_group, Sigma = NULL, jitter = 1e-6) {
  n <- nrow(X_group)
  k <- ncol(X_group)
  if (k < 2) return(X_group)
  if (is.null(Sigma)) {
    Sigma <- stats::cov(X_group)
  }
  Sigma <- .sb_cov_adjust(Sigma, jitter)
  draws <- MASS::mvrnorm(n = n, mu = numeric(k), Sigma = Sigma)
  # centre and normalise each column to unit \ell_2 norm
  draws <- sweep(draws, 2L, colMeans(draws), FUN = "-", check.margin = FALSE)
  norms <- sqrt(colSums(draws^2))
  norms[!is.finite(norms) | norms < sqrt(.Machine$double.eps)] <- 1
  sweep(draws, 2L, norms, "/", check.margin = FALSE)
}

.sb_cache_env <- function(cache) {
  if (is.null(cache)) {
    return(new.env(parent = emptyenv()))
  }
  if (is.environment(cache)) {
    return(cache)
  }
  env <- new.env(parent = emptyenv())
  cache_list <- as.list(cache)
  for (nm in names(cache_list)) {
    assign(nm, cache_list[[nm]], envir = env)
  }
  env
}

.sb_group_diagnostics <- function(original, draws) {
  if (!is.matrix(original)) original <- as.matrix(original)
  size <- ncol(original)
  if (size < 2) {
    return(data.frame(
      mean_abs_corr_orig = NA_real_,
      mean_abs_corr_surrogate = NA_real_,
      mean_abs_corr_cross = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  orig_corr <- stats::cor(original)
  orig_mean <- mean(abs(orig_corr[upper.tri(orig_corr)]))
  if (!is.finite(orig_mean)) orig_mean <- NA_real_
  
  surrogate_means <- vapply(draws, function(mat) {
    cmat <- try(stats::cor(mat), silent = TRUE)
    if (inherits(cmat, "try-error")) return(NA_real_)
    mean(abs(cmat[upper.tri(cmat)]))
  }, numeric(1))
  surrogate_mean <- mean(surrogate_means, na.rm = TRUE)
  if (!is.finite(surrogate_mean)) surrogate_mean <- NA_real_
  
  cross_means <- vapply(draws, function(mat) {
    cmat <- try(stats::cor(original, mat), silent = TRUE)
    if (inherits(cmat, "try-error")) return(NA_real_)
    diag_vals <- diag(cmat)
    mean(abs(diag_vals))
  }, numeric(1))
  cross_mean <- mean(cross_means, na.rm = TRUE)
  if (!is.finite(cross_mean)) cross_mean <- NA_real_
  
  data.frame(
    mean_abs_corr_orig = orig_mean,
    mean_abs_corr_surrogate = surrogate_mean,
    mean_abs_corr_cross = cross_mean,
    stringsAsFactors = FALSE
  )
}

#' Generate correlated design replicates for a set of groups
#'
#' @param X_norm Normalised design matrix.
#' @param groups Correlation structure. Either a list as returned by
#'   [sb_group_variables()] or a vector of group labels matching the columns of
#'   `X_norm`.
#' @param B Number of replicates to generate.
#' @param jitter Numeric value added to covariance diagonals for stability.
#' @param seed Optional integer seed for reproducibility. The seed is scoped via
#'   [withr::with_seed()] so the caller's RNG state is restored on exit.
#' @param use.parallel Logical; when `TRUE`, compute the resampled designs using
#'   a parallel backend when available.
#' @param cache Optional environment or named list used to cache previously
#'   generated surrogates. Passing the same cache across calls reuses draws for
#'   identical groups.
#' @return An object of class `sb_resamples`, i.e. a list of length `B` whose
#'   elements are resampled design matrices. The object exposes per-group
#'   diagnostics in its `"diagnostics"` attribute and returns the cache via the
#'   `"cache"` attribute for reuse.
#' @details When every group has size one (no correlated variables) the function
#'   simply returns `B` copies of `X_norm`. A warning is issued in that situation
#'   so downstream code can avoid mistaking the replicated designs for genuinely
#'   resampled surrogates. The covariance matrices underpinning each correlated
#'   draw are cached in the supplied `cache` environment; reusing the environment
#'   across calls lets `sb_resample_groups()` skip redundant covariance
#'   decompositions for identical groups and speeds up iterative workflows.
#' @export
sb_resample_groups <- function(
    X_norm,
    groups,
    B = 100,
    jitter = 1e-6,
    seed = NULL,
    use.parallel = FALSE,
    cache = NULL
) {
  runner <- function() {
    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
      stop("`B` must be a positive integer.")
    }
    B_int <- as.integer(B)
    if (!isTRUE(all.equal(B, B_int))) {
      warning("Coercing `B` to an integer for resampling.", call. = FALSE)
    }
    X_mat <- if (!is.matrix(X_norm)) as.matrix(X_norm) else X_norm
    if (!is.numeric(X_mat)) {
      storage.mode(X_mat) <- "double"
    }
    original_groups <- groups
    p <- ncol(X_mat)
    as_group_list <- function(g) {
      if (is.list(g)) {
        lapply(g, function(idx) {
          idx <- unique(as.integer(idx))
          idx[!is.na(idx) & idx >= 1 & idx <= p]
        })
      } else {
        if (length(g) != p) {
          stop("`groups` must have length equal to the number of columns when supplied as a vector.")
        }
        f <- as.factor(g)
        split_idx <- split(seq_len(p), f)
        lapply(seq_len(p), function(j) {
          grp <- f[j]
          sort(unname(split_idx[[as.character(grp)]]))
        })
      }
    }
    group_list <- as_group_list(groups)
    if (!length(group_list)) {
      warning("No groups supplied; returning copies of `X_norm`.", call. = FALSE)
    }
    normalised_list <- lapply(group_list, function(idx) sort(unique(idx)))
    has_correlated <- length(normalised_list) && any(vapply(normalised_list, length, integer(1)) >= 2L)
    if (length(group_list) && !has_correlated) {
      warning(
        "All groups are singletons; correlated resampling degenerates to repeated `X_norm`.",
        call. = FALSE
      )
    }
    if (!length(normalised_list)) {
      normalised_list <- list()
    }
    keys <- vapply(normalised_list, function(idx) paste(idx, collapse = ","), character(1))
    unique_pos <- which(keys != "" & !duplicated(keys))
    unique_groups <- normalised_list[unique_pos]

    cache_env <- .sb_cache_env(cache)
    group_draws <- vector("list", length(unique_groups))
    diag_entries <- vector("list", length(unique_groups))
    col_labels <- colnames(X_mat)
    if (is.null(col_labels)) {
      col_labels <- paste0("V", seq_len(ncol(X_mat)))
    }

    if (!length(unique_groups)) {
      resamples <- replicate(B_int, {
        Xb <- X_mat
        attr(Xb, "groups") <- original_groups
        Xb
      }, simplify = FALSE)
      class(resamples) <- c("sb_resamples", "list")
      attr(resamples, "diagnostics") <- data.frame(
        group = character(0),
        size = integer(0),
        regenerated = integer(0),
        cached = logical(0),
        mean_abs_corr_orig = numeric(0),
        mean_abs_corr_surrogate = numeric(0),
        mean_abs_corr_cross = numeric(0)
      )
      attr(resamples, "cache") <- cache_env
      return(resamples)
    }

    for (i in seq_along(unique_groups)) {
      idx <- unique_groups[[i]]
      key <- keys[unique_pos[i]]
      label <- paste(col_labels[idx], collapse = ",")
      if (!nzchar(label)) {
        label <- paste(idx, collapse = ",")
      }

      if (length(idx) < 2) {
        diag_entries[[i]] <- data.frame(
          group = label,
          size = length(idx),
          regenerated = 0L,
          cached = FALSE,
          mean_abs_corr_orig = NA_real_,
          mean_abs_corr_surrogate = NA_real_,
          mean_abs_corr_cross = NA_real_,
          stringsAsFactors = FALSE
        )
        group_draws[i] <- list(NULL)
        next
      }

      stored <- if (exists(key, envir = cache_env, inherits = FALSE)) get(key, envir = cache_env) else NULL
      valid <- !is.null(stored) && is.list(stored$draws) && length(stored$draws) == B_int
      if (isTRUE(valid)) {
        mat0 <- stored$draws[[1]]
        valid <- is.matrix(mat0) && nrow(mat0) == nrow(X_norm) && ncol(mat0) == length(idx)
      }

      if (!valid) {
        Sigma <- stats::cov(X_mat[, idx, drop = FALSE])
        draws <- .sb_parallel_map(
          seq_len(B_int),
          function(b) .sb_sample_group(X_mat[, idx, drop = FALSE], Sigma = Sigma, jitter = jitter),
          use.parallel = use.parallel,
          seed = if (!is.null(seed)) seed + i else NULL
        )
        base_diag <- .sb_group_diagnostics(X_mat[, idx, drop = FALSE], draws)
        stored <- list(draws = draws, diag = base_diag, B = B_int)
        assign(key, stored, envir = cache_env)
        diag_call <- transform(
          base_diag,
          group = label,
          size = length(idx),
          regenerated = B_int,
          cached = FALSE
        )
      } else {
        draws <- stored$draws
        base_diag <- stored$diag
        stored_B <- stored$B
        if (is.null(stored_B) || !is.finite(stored_B)) {
          stored_B <- length(draws)
        }
        diag_call <- transform(
          base_diag,
          group = label,
          size = length(idx),
          regenerated = if (isTRUE(use.parallel)) stored_B else 0L,
          cached = !isTRUE(use.parallel)
        )
      }
      group_draws[[i]] <- draws
      diag_entries[[i]] <- diag_call[, c("group", "size", "regenerated", "cached",
                                         "mean_abs_corr_orig", "mean_abs_corr_surrogate",
                                         "mean_abs_corr_cross")]
    }

    resamples <- vector("list", B_int)
    for (b in seq_len(B_int)) {
      Xb <- X_mat
      for (i in seq_along(unique_groups)) {
        idx <- unique_groups[[i]]
        draws <- group_draws[[i]]
        if (is.null(draws)) next
        Xb[, idx] <- draws[[b]]
      }
      attr(Xb, "groups") <- original_groups
      resamples[[b]] <- Xb
    }

    class(resamples) <- c("sb_resamples", "list")
    diag_df <- if (length(diag_entries)) do.call(rbind, diag_entries) else data.frame(
      group = character(0),
      size = integer(0),
      regenerated = integer(0),
      cached = logical(0),
      mean_abs_corr_orig = numeric(0),
      mean_abs_corr_surrogate = numeric(0),
      mean_abs_corr_cross = numeric(0)
    )
    attr(resamples, "diagnostics") <- diag_df
    attr(resamples, "cache") <- cache_env
    resamples
  }

  if (is.null(seed)) {
    runner()
  } else {
    withr::with_seed(seed, runner())
  }
}

#' Apply a selector to a collection of resampled designs
#'
#' @param X_norm Normalised design matrix.
#' @param resamples List of matrices returned by [sb_resample_groups()].
#' @param Y Numeric response.
#' @param selector Variable-selection routine; function or character string. If
#'   it is a function, the selector name should be added as the fun.name attribute.
#' @param ... Extra arguments passed to the selector.
#' @return A numeric matrix of coefficients with one column per resample.
#' @param keep_template Logical; when `TRUE` (default) the first column stores the
#'   coefficients fitted on `X_norm` before any resampling. This avoids
#'   recomputing the selector on the original data and keeps the baseline fit
#'   available for diagnostics.
#' @export
sb_apply_selector_manual <- function(X_norm, resamples, Y, selector, ..., keep_template = TRUE) {
  fun <- if (is.character(selector)) get(selector, mode = "function", envir = parent.frame()) else selector
  template <- fun(X_norm, Y, ...)
  if (!is.numeric(template)) {
    stop("`selector` must return a numeric vector of coefficients.")
  }
  template_names <- names(template)
  if (is.null(template_names)) {
    stop("`selector` must return named coefficients so results can be aligned.")
  }
  coef_len <- length(template)
  n_resamples <- length(resamples)
  offset <- if (isTRUE(keep_template)) 1L else 0L
  total_cols <- n_resamples + offset
  if (total_cols == 0L) {
    mat <- matrix(template, ncol = 1L)
    rownames(mat) <- template_names
    colnames(mat) <- "sim0"
    return(mat)
  }
  out <- matrix(NA_real_, coef_len, total_cols)
  base_colnames <- paste0("sim", seq_len(n_resamples))
  if (isTRUE(keep_template)) {
    colnames(out) <- c("sim0", base_colnames)
    out[, 1L] <- template
  } else {
    colnames(out) <- base_colnames
  }
  rownames(out) <- template_names

  align_output <- function(beta, replicate_id) {
    if (!is.numeric(beta)) {
      stop(sprintf("Selector returned non-numeric coefficients at replicate %d.", replicate_id))
    }
    nm <- names(beta)
    if (is.null(nm)) {
      stop(sprintf("Selector returned unnamed coefficients at replicate %d.", replicate_id))
    }
    idx <- match(template_names, nm)
    if (anyNA(idx)) {
      missing <- template_names[is.na(idx)]
      stop(sprintf(
        "Selector output is missing coefficients for: %s (replicate %d).",
        paste(missing, collapse = ", "),
        replicate_id
      ))
    }
    as.numeric(beta[idx])
  }

  for (b in seq_len(n_resamples)) {
    target_col <- b + offset
    if (isTRUE(keep_template) && b == 1L && isTRUE(all.equal(resamples[[b]], X_norm))) {
      out[, target_col] <- template
      next
    }
    beta <- fun(resamples[[b]], Y, ...)
    out[, target_col] <- align_output(beta, b)
  }
  out
}

#' Compute selection frequencies from coefficient paths
#'
#' @param coef_matrix Matrix produced by [sb_apply_selector_manual()].
#' @param version Either `"glmnet"` (first row is intercept) or `"lars"`.
#' @param threshold Coefficients with absolute value below this threshold are
#'   treated as zero.
#' @return Numeric vector of selection frequencies.
#' @export
sb_selection_frequency <- function(coef_matrix, version = c("glmnet", "lars"), threshold = 1e-4) {
  version <- match.arg(version)
  mat <- coef_matrix
  if (version == "glmnet") {
    mat <- mat[-1, , drop = FALSE]
  }
  if (!is.matrix(mat)) {
    mat <- matrix(mat, ncol = 1L)
  }
  freq <- rowMeans(abs(mat) > threshold)
  names(freq) <- rownames(mat)
  freq
}

# Internal: derive c0 sequence from correlation values
.sb_c0_sequence <- function(corr_mat, step.num = 0.1, steps.seq = NULL) {
  vals <- abs(corr_mat[upper.tri(corr_mat)])
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(numeric(0))
  vals <- sort(unique(vals))
  if (!is.null(steps.seq)) {
    seq <- sort(unique(pmax(pmin(steps.seq, 1), 0)), decreasing = TRUE)
    seq[seq > 0 & seq < 1]
  } else {
    step <- max(step.num, if (length(vals) > 2) 1 / (length(vals) - 2) else step.num)
    probs <- rev(unique(c(0, seq(from = 0, to = 1, by = step), 1)))
    probs <- probs[probs > 0 & probs < 1]
    stats::quantile(vals, probs, names = FALSE, type = 1)
  }
}

#' SelectBoost for beta-regression models
#'
#' @description
#' `sb_beta()` orchestrates all SelectBoost stages—normalisation, correlation
#' analysis, grouping, correlated resampling, and stability tallying—while using
#' the beta-regression selectors provided by this package. It can operate on
#' point-valued or interval-valued responses and automatically squeezes the
#' outcome into `(0, 1)` unless instructed otherwise.
#'
#' @inheritParams sb_resample_groups
#' @param X Numeric design matrix. Coerced with [as.matrix()] and normalised via
#'   [sb_normalize()].
#' @param Y Numeric response vector. Values are squeezed to the open unit
#'   interval with the standard SelectBoost transformation unless `squeeze =
#'   FALSE`. Optional when interval bounds are supplied.
#' @param selector Selection routine. Defaults to [betareg_step_aic()]. 
#'   Function or character string. If it is a function, the selector name 
#'   should be added as the fun.name attribute.
#' @param corrfunc Correlation function passed to [sb_compute_corr()].
#' @param step.num Step length for the automatically generated `c0` grid.
#' @param steps.seq Optional user-supplied grid of absolute correlation
#'   thresholds.
#' @param version Either `"glmnet"` (intercept in first row) or `"lars"`.
#' @param squeeze Logical; ensure the response lies in `(0, 1)`.
#' @param use.parallel Logical; enable parallel resampling and selector fits
#'   when supported by the current R session.
#' @param verbose Logical; emit progress messages.
#' @param threshold Numeric tolerance for considering a coefficient selected.
#' @param interval Interval-resampling mode: `"none"` reuses `Y`, whereas
#'   `"uniform"` and `"midpoint"` draw pseudo-responses between `Y_low` and
#'   `Y_high` for each replicate.
#' @param Y_low,Y_high Interval bounds in `[0, 1]` paired with the rows of `X`
#'   when `interval` is not `"none"`.
#' @param ... Additional arguments forwarded to `selector`.
#' @details
#' The returned object carries a rich set of attributes:
#'
#' * `"c0.seq"` – the grid of absolute-correlation thresholds explored during
#'   resampling.
#' * `"steps.seq"` – the raw sequence (if any) used to construct the grid.
#' * `"selector"` – the selector identifier (function name or expression).
#' * `"B"` – number of resampled designs passed to the selector.
#' * `"interval"` – the interval sampling mode (`"none"`, `"uniform"`, or
#'   `"midpoint"`).
#' * `"resample_diagnostics"` – per-threshold data frames with summary
#'   statistics on the cached correlated draws.
#'
#' These attributes mirror the historical SelectBoost beta implementation so the
#' object can be consumed by existing plotting and reporting utilities.
#'
#' @return Matrix of selection frequencies with one row per `c0` level and class
#'   `"sb_beta"`. See *Details* for the recorded attributes.
#' @examples
#' set.seed(42)
#' sim <- simulation_DATA.beta(n = 80, p = 4, s = 2)
#' # increase B for real applications
#' res <- sb_beta(sim$X, sim$Y, B = 5)
#' res
#' @export
sb_beta <- function(
    X,
    Y = NULL,
    selector = betareg_step_aic,
    corrfunc = "cor",
    B = 100,
    step.num = 0.1,
    steps.seq = NULL,
    version = c("glmnet", "lars"),
    squeeze = TRUE,
    use.parallel = FALSE,
    seed = NULL,
    verbose = FALSE,
    threshold = 1e-4,
    interval = c("none", "uniform", "midpoint"),
    Y_low = NULL,
    Y_high = NULL,
    ...
) {
  version_local  <- match.arg(version,  c("glmnet", "lars"))
  interval_local <- match.arg(interval, c("none", "uniform", "midpoint"))

  runner <- function() {
    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
      stop("`B` must be a positive integer.")
    }
    B_int <- as.integer(B)
    if (!isTRUE(all.equal(B, B_int))) {
      warning("Coercing `B` to an integer for sb_beta().", call. = FALSE)
    }
    X_norm <- sb_normalize(X)
    n <- nrow(X_norm)
    if (!is.null(Y) && length(Y) != n) {
      stop("`Y` must have length equal to nrow(X).")
    }
    if (interval_local == "none" && is.null(Y)) {
      stop("`Y` must be supplied when `interval = \"none\"`.")
    }
    y_input <- if (!is.null(Y)) as.numeric(Y) else rep(NA_real_, n)
    if (!all(is.na(y_input))) {
      if (any(!is.finite(y_input[!is.na(y_input)]))) {
        stop("`Y` must contain only finite values.")
      }
    }

    clamp01 <- function(v) {
      pmin(pmax(v, 0), 1)
    }

    if (interval_local != "none") {
      if (is.null(Y_low) || is.null(Y_high)) {
        stop("`Y_low` and `Y_high` must be supplied when `interval` is not `\"none\"`.")
      }
      if (length(Y_low) != n || length(Y_high) != n) {
        stop("`Y_low` and `Y_high` must have length equal to nrow(X).")
      }
      y_low <- clamp01(as.numeric(Y_low))
      y_high <- clamp01(as.numeric(Y_high))
      if (any(!is.finite(y_low) | !is.finite(y_high))) {
        stop("`Y_low` and `Y_high` must contain finite values.")
      }
      if (any(y_low > y_high)) {
        stop("Each element of `Y_low` must be <= the corresponding element of `Y_high`.")
      }
      y_base <- if (!is.null(Y)) y_input else 0.5 * (y_low + y_high)
      y_base <- clamp01(y_base)
    } else {
      y_base <- y_input
    }

    if (any(is.na(y_base))) {
      stop("`Y` contains missing values.")
    }

    if (squeeze) {
      y_base <- .sqz01(y_base)
    } else if (any(y_base <= 0 | y_base >= 1)) {
      stop("`Y` must lie in (0, 1) when `squeeze = FALSE`.")
    }

    sample_interval <- function() {
      if (interval_local == "midpoint") {
        vals <- 0.5 * (y_low + y_high)
      } else {
        u <- stats::runif(n)
        vals <- y_low + u * (y_high - y_low)
      }
      vals <- clamp01(vals)
      if (squeeze) {
        .sqz01(vals)
      } else {
        pmin(pmax(vals, .Machine$double.eps), 1 - .Machine$double.eps)
      }
    }

    B_rep <- max(1L, B_int)
    y_draws <- if (interval_local == "none") {
      rep(list(y_base), B_rep)
    } else {
      replicate(B_rep, sample_interval(), simplify = FALSE)
    }
    y_template <- y_base
    corr_mat <- sb_compute_corr(X_norm, corrfunc = corrfunc)
    c0_inner <- .sb_c0_sequence(corr_mat, step.num = step.num, steps.seq = steps.seq)
    c0_levels <- unique(c(1, c0_inner, 0))
    res_list <- vector("list", length(c0_levels))
    names(res_list) <- sprintf("c0 = %.3f", c0_levels)
    diag_list <- vector("list", length(c0_levels))
    if (verbose) {
      message("Running sb_beta over ", length(c0_levels), " correlation thresholds")
    }
    fun <- if (is.character(selector)) get(selector, mode = "function", envir = parent.frame()) else selector
    template <- fun(X_norm, y_template, ...)
    if (!is.numeric(template)) {
      stop("`selector` must return numeric coefficients.")
    }
    coef_names <- names(template)
    if (is.null(coef_names)) {
      stop("`selector` must return named coefficients for alignment.")
    }
    empty_diag <- data.frame(
      group = character(0),
      size = integer(0),
      regenerated = integer(0),
      cached = logical(0),
      mean_abs_corr_orig = numeric(0),
      mean_abs_corr_surrogate = numeric(0),
      mean_abs_corr_cross = numeric(0)
    )
    resample_cache <- new.env(parent = emptyenv())

    align_selector <- function(beta, replicate_id) {
      if (!is.numeric(beta)) {
        stop(sprintf("Selector returned non-numeric coefficients at replicate %d.", replicate_id))
      }
      nm <- names(beta)
      if (is.null(nm)) {
        stop(sprintf("Selector returned unnamed coefficients at replicate %d.", replicate_id))
      }
      idx <- match(coef_names, nm)
      if (anyNA(idx)) {
        missing <- coef_names[is.na(idx)]
        stop(sprintf(
          "Selector output is missing coefficients for: %s (replicate %d).",
          paste(missing, collapse = ", "),
          replicate_id
        ))
      }
      as.numeric(beta[idx])
    }

    for (i in seq_along(c0_levels)) {
      c0 <- c0_levels[i]
      groups <- sb_group_variables(corr_mat, c0)
      if (verbose) message(sprintf("c0 = %.3f", c0))
      big_groups <- Filter(function(idx) length(idx) >= 2, groups)
      if (!length(big_groups)) {
        diag_list[[i]] <- empty_diag
        if (interval_local == "none") {
          coef_mat <- matrix(template, ncol = 1L)
          rownames(coef_mat) <- coef_names
          colnames(coef_mat) <- "sim1"
        } else {
          resamples <- replicate(B_rep, {
            mat <- X_norm
            attr(mat, "groups") <- groups
            mat
          }, simplify = FALSE)
          y_vals <- y_draws
          coef_vals <- .sb_parallel_map(
            seq_along(resamples),
            function(b) align_selector(fun(resamples[[b]], y_vals[[b]], ...), b),
            use.parallel = use.parallel,
            seed = if (!is.null(seed)) seed + 1000L + i
          )
          coef_mat <- do.call(cbind, coef_vals)
          rownames(coef_mat) <- coef_names
          colnames(coef_mat) <- paste0("sim", seq_along(resamples))
        }
      } else {
        resample_input <- if (length(big_groups)) big_groups else list()
        resamples <- sb_resample_groups(
          X_norm,
          resample_input,
          B = B_int,
          jitter = 1e-6,
          seed = if (!is.null(seed)) seed + i,
          use.parallel = use.parallel,
          cache = resample_cache
        )
        resample_cache <- attr(resamples, "cache")
        diag_list[[i]] <- attr(resamples, "diagnostics")
        y_vals <- if (interval_local == "none") {
          rep(list(y_template), length(resamples))
        } else {
          y_draws
        }
        coef_vals <- .sb_parallel_map(
          seq_along(resamples),
          function(b) align_selector(fun(resamples[[b]], y_vals[[b]], ...), b),
          use.parallel = use.parallel,
          seed = if (!is.null(seed)) seed + 1000L + i
        )
        coef_mat <- do.call(cbind, coef_vals)
        rownames(coef_mat) <- coef_names
        colnames(coef_mat) <- paste0("sim", seq_along(resamples))
      }
      freq <- sb_selection_frequency(coef_mat, version = version_local, threshold = threshold)
      res_list[[i]] <- freq
      if (is.null(diag_list[[i]])) {
        diag_list[[i]] <- empty_diag
      }
    }
    names(diag_list) <- names(res_list)
    freq_mat <- do.call(rbind, res_list)
    attr(freq_mat, "c0.seq") <- c0_levels
    attr(freq_mat, "steps.seq") <- c0_inner
    attr(freq_mat, "B") <- B_int
    attr(freq_mat, "selector") <- if (is.character(selector)) selector else {
      if(is.null(attr(selector,"fun.name"))){
        sel_expr <- substitute(selector)
        if (is.symbol(sel_expr)) as.character(sel_expr) else paste(deparse(sel_expr), collapse = "")
      } else {
        attr(selector,"fun.name")
      }
    }
    attr(freq_mat, "resample_diagnostics") <- diag_list
    attr(freq_mat, "interval") <- interval_local
    class(freq_mat) <- c("sb_beta", class(freq_mat))
    freq_mat
  }

  if (is.null(seed)) {
    runner()
  } else {
    withr::with_seed(seed, runner())
  }
}
