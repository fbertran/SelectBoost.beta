#' Run all selectors once on a dataset
#'
#' Convenience wrapper that runs AIC/BIC/AICc stepwise, GAMLSS LASSO (and ENet
#' when available), and the pure glmnet IRLS selector, then collates coefficients
#' into a long table for comparison. Observations containing `NA` in either `X`
#' or `Y` are removed prior to fitting. Column names are temporarily shortened
#' to satisfy selector requirements and avoid clashes; the outputs remap them to
#' the original labels before returning so the reported variables always match
#' the input design.
#'
#' @inheritParams betareg_step_aic
#' @param include_enet Logical; include ENet if `gamlss.lasso` is installed.
#'
#' @return A list with:
#' \describe{
#'   \item{coefs}{Named coefficient vectors for each selector.}
#'   \item{table}{Long data frame with columns `selector`, `variable`, `coef`, `selected`.}
#' }
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[, 1])
#' Y <- rbeta(100, Y * 30, (1 - Y) * 30)
#' single <- compare_selectors_single(X, Y, include_enet = FALSE)
#' head(single$table)
#' @export
compare_selectors_single <- function(X, Y, include_enet = TRUE) {
  Y <- as.numeric(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (length(Y) != nrow(X)) {
    stop("`Y` must have length equal to `nrow(X)`.")
  }
  
  keep <- stats::complete.cases(X) & !is.na(Y)
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    Y <- Y[keep]
  }
  if (!nrow(X)) {
    stop("No complete cases available after removing missing values.")
  }
  
  ss <- .shorten_colnames(X)
  Xs <- ss$X; map <- ss$map
  
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  sels <- list(
    AIC   = betareg_step_aic,
    BIC   = betareg_step_bic,
    AICc  = betareg_step_aicc
  )
  has_gamlss <- requireNamespace("gamlss", quietly = TRUE)
  if (has_gamlss) {
    sels$LASSO <- betareg_lasso_gamlss
  }
  sels$GLMNET <- function(X, Y) betareg_glmnet(X, Y, choose = "bic", n_iter = 5, prestandardize = TRUE)
  if (has_gamlss && include_enet && requireNamespace("gamlss.lasso", quietly = TRUE)) {
    sels$ENET <- function(X, Y) betareg_enet_gamlss(X, Y, method = "IC", ICpen = "BIC", alpha = 0.8)
  }

#  sels$AIC(Xs,Y)
#  sels$BIC(Xs,Y)
#  sels$AICc(Xs,Y)
#  sels$LASSO(Xs,Y)
#  sels$GLMNET(Xs,Y)
#  sels$ENET(Xs,Y)
  
  coefs <- lapply(sels, function(f) suppressWarnings(f(Xs, Y)))
  # remap names back to original for outputs
  coefs <- lapply(coefs, function(b) {
    out <- b
    nm <- names(b)
    # keep intercept, remap others via 'map'
    if (length(nm) > 1) {
      short <- nm[-1]
      orig  <- unname(map[short])
      names(out) <- c("(Intercept)", orig)
    }
    out
  })
  vars <- unname(map[colnames(Xs)])
  tab <- do.call(rbind, lapply(names(coefs), function(nm) {
    b <- coefs[[nm]][vars]
    data.frame(selector = nm, variable = vars, coef = as.numeric(b), selected = b != 0)
  }))
  list(coefs = coefs, table = tab)
}

#' Bootstrap selection frequencies across selectors
#'
#' Bootstraps the dataset `B` times and records how often each variable is
#' selected by each selector. Observations containing `NA` in either `X` or `Y`
#' are removed prior to resampling. Column names are abbreviated internally and
#' mapped back to the originals in the output just like in
#' [compare_selectors_single()].
#'
#' @inheritParams compare_selectors_single
#' @param B Number of bootstrap replications.
#' @param seed Optional RNG seed.
#'
#' @return Long data frame with columns `selector`, `variable`, `freq` in `[0,1]`,
#'   `n_success`, and `n_fail`. The `freq` column reports the share of bootstrap
#'   replicates where a variable was selected by the corresponding selector.
#'   Values near 1 signal high stability whereas small values indicate weak
#'   evidence. `n_success` counts the successful fits contributing to the
#'   frequency estimate (excluding failed replicates), while `n_fail` records the
#'   number of unsuccessful fits. A `"failures"` attribute attached to the
#'   returned data frame lists the replicate indices and messages for any
#'   encountered errors.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[, 1])
#' Y <- rbeta(100, Y * 30, (1 - Y) * 30)
#' freq <- compare_selectors_bootstrap(X, Y, B = 10, include_enet = FALSE)
#' head(freq)
#' subset(freq, freq > 0.8)
#'
#' # Increase `B` until the reported frequencies stabilise. For example,
#' # freq_big <- compare_selectors_bootstrap(X, Y, B = 200, include_enet = FALSE)
#' # stats::aggregate(freq ~ selector, freq_big, summary)
#'
#' @export
compare_selectors_bootstrap <- function(X, Y, B = 50, include_enet = TRUE, seed = NULL) {
  runner <- function() {
    Y <- as.numeric(Y)
    if (!is.matrix(X)) X <- as.matrix(X)
    if (length(Y) != nrow(X)) {
      stop("`Y` must have length equal to `nrow(X)`.")
    }
    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1) {
      stop("`B` must be a positive integer.")
    }
    B_int <- as.integer(B)
    if (!isTRUE(all.equal(B, B_int))) {
      warning("Coercing `B` to an integer for bootstrapping.", call. = FALSE)
    }

    keep <- stats::complete.cases(X) & !is.na(Y)
    if (!all(keep)) {
      X <- X[keep, , drop = FALSE]
      Y <- Y[keep]
    }
    if (!nrow(X)) {
      stop("No complete cases available after removing missing values.")
    }

    if (!is.matrix(X)) X <- as.matrix(X)
    if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

    ss <- .shorten_colnames(X)
    Xs <- ss$X
    map <- ss$map

    n <- nrow(Xs)
    short_vars <- colnames(Xs)
    vars <- unname(map[short_vars])

    base <- list(
      AIC   = betareg_step_aic,
      BIC   = betareg_step_bic,
      AICc  = betareg_step_aicc
    )
    has_gamlss <- requireNamespace("gamlss", quietly = TRUE)
    if (has_gamlss) {
      base$LASSO <- betareg_lasso_gamlss
    }
    base$GLMNET <- function(X, Y) betareg_glmnet(X, Y, choose = "bic", n_iter = 4, prestandardize = TRUE)
    if (has_gamlss && include_enet && requireNamespace("gamlss.lasso", quietly = TRUE)) {
      base$ENET <- function(X, Y) betareg_enet_gamlss(X, Y, method = "IC", ICpen = "BIC", alpha = 0.8)
    }

    align_coefficients <- function(beta, selector_name, replicate_id) {
      if (!is.numeric(beta)) {
        stop(sprintf("Selector `%s` returned a non-numeric coefficient vector at replicate %d.", selector_name, replicate_id))
      }
      if (is.null(names(beta))) {
        stop(sprintf("Selector `%s` returned unnamed coefficients at replicate %d.", selector_name, replicate_id))
      }
      nm <- names(beta)
      if (length(nm) > 1L) {
        mapped <- unname(map[nm[!(grepl("(Intercept)", nm) | grepl("^phi\\|", nm))]])
        if (any(is.na(mapped))) {
          stop(sprintf(
            "Selector `%s` produced coefficients that could not be mapped back to the original names at replicate %d.",
            selector_name,
            replicate_id
          ))
        }
        nm[!(grepl("(Intercept)", nm) | grepl("^phi\\|", nm))] <- mapped
      }
      names(beta) <- nm
      idx <- match(vars, names(beta))
      if (anyNA(idx)) {
        missing <- vars[is.na(idx)]
        stop(sprintf(
          "Selector `%s` is missing coefficients for: %s (replicate %d).",
          selector_name,
          paste(missing, collapse = ", "),
          replicate_id
        ))
      }
      beta[idx]
    }

    failures <- setNames(vector("list", length(base)), names(base))

    freq_list <- lapply(names(base), function(nm) {
      f <- base[[nm]]
      sel <- matrix(NA_real_, B_int, length(vars))
      colnames(sel) <- vars
      fail_log <- list()
      for (b in seq_len(B_int)) {
        idx <- sample.int(n, n, replace = TRUE)
        beta <- tryCatch(
          f(Xs[idx, , drop = FALSE], Y[idx]),
          error = function(e) {
            structure(list(condition = e), class = "sb_selector_failure")
          }
        )
        if (inherits(beta, "sb_selector_failure")) {
          fail_log[[length(fail_log) + 1L]] <- data.frame(
            replicate = b,
            message = conditionMessage(beta$condition),
            stringsAsFactors = FALSE
          )
          next
        }

        aligned <- tryCatch(
          align_coefficients(beta, nm, b),
          error = function(e) {
            structure(list(condition = e), class = "sb_selector_failure")
          }
        )
        if (inherits(aligned, "sb_selector_failure")) {
          fail_log[[length(fail_log) + 1L]] <- data.frame(
            replicate = b,
            message = conditionMessage(aligned$condition),
            stringsAsFactors = FALSE
          )
          next
        }

        sel[b, ] <- as.numeric(aligned != 0)
      }
      successes <- as.integer(colSums(!is.na(sel)))
      freq <- colMeans(sel, na.rm = TRUE)
      freq[successes == 0] <- NA_real_
      failures[[nm]] <<- if (length(fail_log)) {
        do.call(rbind, fail_log)
      } else {
        data.frame(replicate = integer(0), message = character(0))
      }
      data.frame(
        selector = nm,
        variable = vars,
        freq = freq,
        n_success = successes,
        n_fail = B_int - successes
      )
    })

    result <- do.call(rbind, freq_list)
    attr(result, "failures") <- failures
    result
  }

  if (is.null(seed)) {
    runner()
  } else {
    withr::with_seed(seed, runner())
  }
}


#' Merge single-run results and bootstrap frequencies
#'
#' @param single_tab Data frame returned in `compare_selectors_single()[["table"]]`.
#' @param freq_tab Optional frequency table from [compare_selectors_bootstrap()].
#' @return Merged data frame.
#' @examples
#' single_tab <- data.frame(
#'   selector = rep(c("AIC", "BIC"), each = 3),
#'   variable = rep(paste0("x", 1:3), times = 2),
#'   coef = c(0.5, 0, -0.2, 0.6, 0.1, -0.3)
#' )
#' single_tab$selected <- single_tab$coef != 0
#' freq_tab <- data.frame(
#'   selector = rep(c("AIC", "BIC"), each = 3),
#'   variable = rep(paste0("x", 1:3), times = 2),
#'   freq = c(0.9, 0.15, 0.4, 0.85, 0.3, 0.25)
#' )
#' compare_table(single_tab, freq_tab)
#' @export
compare_table <- function(single_tab, freq_tab = NULL) {
  if (!is.null(freq_tab)) {
    merge(single_tab, freq_tab, by = c("selector","variable"), all = TRUE)
  } else single_tab
}


#' Side-by-side coefficient heatmap
#'
#' Visual comparison of coefficients returned by each selector. Requires `ggplot2`.
#'
#' @param single_tab Data frame as returned by `compare_selectors_single()[["table"]]`.
#' @return A `ggplot` object when `ggplot2` is available; otherwise draws a base R image.
#' @importFrom rlang .data
#' @examples
#' demo_tab <- data.frame(
#'   selector = rep(c("AIC", "BIC"), each = 3),
#'   variable = rep(paste0("x", 1:3), times = 2),
#'   coef = c(0.6, 0, -0.2, 0.55, 0.05, -0.3)
#' )
#' demo_tab$selected <- demo_tab$coef != 0
#' plot_compare_coeff(demo_tab)
#' @export
plot_compare_coeff <- function(single_tab) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(
      single_tab,
      ggplot2::aes(x = .data$variable, y = .data$selector, fill = .data$coef)
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(.data$selected, sprintf("%.2f", .data$coef), "")),
        size = 3
      ) +
      ggplot2::scale_fill_gradient2(name = "coef", na.value = "grey90") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    tab <- reshape(single_tab[, c("selector","variable","coef")], timevar = "variable",
                   idvar = "selector", direction = "wide")
    mat <- as.matrix(tab[,-1, drop=FALSE]); rownames(mat) <- tab[,1]
    image(t(mat[nrow(mat):1, , drop=FALSE]), axes=FALSE, main="Coefficients (no ggplot2)")
    box()
  }
}

#' Side-by-side selection-frequency heatmap
#'
#' Visual comparison of bootstrap selection frequencies by selector. Requires `ggplot2`.
#'
#' @param freq_tab Data frame as returned by [compare_selectors_bootstrap()].
#' @return A `ggplot` object when `ggplot2` is available; otherwise draws a base R image.
#' @importFrom rlang .data
#' @examples
#' freq_tab <- data.frame(
#'   selector = rep(c("AIC", "BIC"), each = 3),
#'   variable = rep(paste0("x", 1:3), times = 2),
#'   freq = c(0.85, 0.2, 0.45, 0.75, 0.35, 0.3)
#' )
#' plot_compare_freq(freq_tab)
#' @export
plot_compare_freq <- function(freq_tab) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(
      freq_tab,
      ggplot2::aes(x = .data$variable, y = .data$selector, fill = .data$freq)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(limits = c(0,1), name = "freq") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    tab <- reshape(freq_tab[, c("selector","variable","freq")], timevar = "variable",
                   idvar = "selector", direction = "wide")
    mat <- as.matrix(tab[,-1, drop=FALSE]); rownames(mat) <- tab[,1]
    image(t(mat[nrow(mat):1, , drop=FALSE]), axes=FALSE, main="Selection freq (no ggplot2)")
    box()
  }
}
