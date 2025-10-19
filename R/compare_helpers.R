#' Run all selectors once on a dataset
#'
#' Convenience wrapper that runs AIC/BIC/AICc stepwise, GAMLSS LASSO (and ENet
#' when available), and the pure glmnet IRLS selector, then collates coefficients
#' into a long table for comparison.
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
#' \dontrun{
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1])
#' Y <- rbeta(100, Y*30, (1-Y)*30); out <- compare_selectors_single(X, Y)
#' head(out$table)
#' }
#' @export
compare_selectors_single <- function(X, Y, include_enet = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  sels <- list(
    AIC   = betareg_step_aic,
    BIC   = betareg_step_bic,
    AICc  = betareg_step_aicc,
    LASSO = betareg_lasso_gamlss,
    GLMNET= function(X,Y) betareg_glmnet(X,Y, choose="bic", n_iter=5, prestandardize=TRUE)
  )
  if (include_enet && requireNamespace("gamlss.lasso", quietly = TRUE)) {
    sels$ENET <- function(X,Y) betareg_enet_gamlss(X, Y, method="IC", ICpen="BIC", alpha=0.8)
  }
  coefs <- lapply(sels, function(f) f(X,Y))
  vars <- colnames(X)
  tab <- do.call(rbind, lapply(names(coefs), function(nm) {
    b <- coefs[[nm]][vars]
    data.frame(selector = nm, variable = vars, coef = as.numeric(b), selected = b != 0)
  }))
  list(coefs = coefs, table = tab)
}

#' Bootstrap selection frequencies across selectors
#'
#' Bootstraps the dataset `B` times and records how often each variable is
#' selected by each selector.
#'
#' @inheritParams compare_selectors_single
#' @param B Number of bootstrap replications.
#' @param seed Optional RNG seed.
#'
#' @return Long data frame with columns `selector`, `variable`, `freq` in `[0,1]`.
#' @examples
#' \dontrun{
#' set.seed(1); X <- matrix(rnorm(300), 100, 3); Y <- plogis(X[,1])
#' Y <- rbeta(100, Y*30, (1-Y)*30); fr <- compare_selectors_bootstrap(X, Y, B = 20)
#' head(fr)
#' }
#' @export
compare_selectors_bootstrap <- function(X, Y, B = 50, include_enet = TRUE, seed = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X); vars <- colnames(X)
  base <- list(
    AIC   = betareg_step_aic,
    BIC   = betareg_step_bic,
    AICc  = betareg_step_aicc,
    LASSO = betareg_lasso_gamlss,
    GLMNET= function(X,Y) betareg_glmnet(X,Y, choose="bic", n_iter=4, prestandardize=TRUE)
  )
  if (include_enet && requireNamespace("gamlss.lasso", quietly = TRUE)) {
    base$ENET <- function(X,Y) betareg_enet_gamlss(X, Y, method="IC", ICpen="BIC", alpha=0.8)
  }
  freq_list <- lapply(names(base), function(nm) {
    f <- base[[nm]]
    sel <- matrix(FALSE, B, length(vars)); colnames(sel) <- vars
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      beta <- f(X[idx,,drop=FALSE], Y[idx])
      sel[b,] <- beta[vars] != 0
    }
    data.frame(selector = nm, variable = vars, freq = colMeans(sel))
  })
  do.call(rbind, freq_list)
}


#' Merge single-run results and bootstrap frequencies
#'
#' @param single_tab Data frame returned in `compare_selectors_single()[["table"]]`.
#' @param freq_tab Optional frequency table from [compare_selectors_bootstrap()].
#' @return Merged data frame.
#' @examples
#' \dontrun{
#' # tbl <- compare_table(single$table, freq)
#' }
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
#' \dontrun{
#' # plot_compare_coeff(single$table)
#' }
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
#' \dontrun{
#' # plot_compare_freq(freq)
#' }
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
