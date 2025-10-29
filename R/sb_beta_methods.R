#' User-friendly methods for `sb_beta()` results
#'
#' These S3 helpers make it easier to inspect and visualise the
#' correlation-threshold grid returned by [sb_beta()]. They surface the stored
#' attributes, reshape the selection frequencies into tidy summaries, and produce
#' quick ggplot2 visualisations for interactive use.
#'
#' @rdname sb_beta_methods
#' @name sb_beta_methods
#' @aliases print.sb_beta, summary.sb_beta, print.summary.sb_beta, autoplot.sb_beta, NULL
#'
#' @param x,object An object of class `sb_beta`.
#' @param digits Number of decimal places to display when printing.
#' @param n Number of rows to show from the summary table when printing.
#' @param variables Optional character vector of variables to retain in the
#'   plotted output.
#' @param ... Additional arguments passed on to lower-level methods.
#' @return `summary.sb_beta()` returns an object of class `summary.sb_beta`
#'   containing a tidy data frame of selection frequencies. The plotting and
#'   printing methods are invoked for their side effects and return the input
#'   object invisibly.
#' @examples
#' set.seed(42)
#' sim <- simulation_DATA.beta(n = 50, p = 4, s = 2)
#' fit <- sb_beta(sim$X, sim$Y, B = 5, step.num = 0.5)
#' print(fit)
#' summary(fit)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   autoplot.sb_beta(fit)
#' }
#' 
NULL

#' Print of `sb_beta()` results
#' @rdname sb_beta_methods
#' @export
print.sb_beta <- function(x, digits = 3, ...) {
  freq_mat <- as.matrix(x)
  freq_mat <- unclass(freq_mat)
  selector <- attr(x, "selector", exact = TRUE)
  B <- attr(x, "B", exact = TRUE)
  c0_levels <- attr(x, "c0.seq", exact = TRUE)
  steps_seq <- attr(x, "steps.seq", exact = TRUE)
  interval <- attr(x, "interval", exact = TRUE)
  cat("SelectBoost beta selection frequencies\n")
  if (!is.null(selector)) {
    cat("Selector: ", selector, "\n", sep = "")
  }
  if (!is.null(B)) {
    cat("Resamples per threshold: ", B, "\n", sep = "")
  }
  if (!is.null(interval)) {
    cat("Interval mode: ", interval, "\n", sep = "")
  }
  if (!is.null(c0_levels)) {
    cat("c0 grid: ", paste(format(round(c0_levels, digits), trim = TRUE), collapse = ", "), "\n", sep = "")
  }
  if (!is.null(steps_seq) && length(steps_seq)) {
    cat("Inner thresholds: ", paste(format(round(steps_seq, digits), trim = TRUE), collapse = ", "), "\n", sep = "")
  }
  if (length(freq_mat)) {
    matrix_to_print <- round(freq_mat, digits)
    print(matrix_to_print)
  } else {
    cat("(no frequencies stored)\n")
  }
  invisible(x)
}

#' Summary of `sb_beta()` results
#' @rdname sb_beta_methods
#' @export
summary.sb_beta <- function(object, ...) {
  freq_mat <- as.matrix(object)
  selector <- attr(object, "selector", exact = TRUE)
  B <- attr(object, "B", exact = TRUE)
  c0_levels <- attr(object, "c0.seq", exact = TRUE)
  interval <- attr(object, "interval", exact = TRUE)
  steps_seq <- attr(object, "steps.seq", exact = TRUE)
  if (is.null(c0_levels)) {
    c0_levels <- seq_len(nrow(freq_mat))
  }
  var_names <- colnames(freq_mat)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(ncol(freq_mat)))
  }
  summary_df <- data.frame(
    c0 = rep(c0_levels, each = ncol(freq_mat)),
    variable = rep(var_names, times = nrow(freq_mat)),
    frequency = as.numeric(freq_mat),
    row.names = NULL
  )
  res <- list(
    selector = selector,
    B = B,
    c0 = c0_levels,
    steps.seq = steps_seq,
    interval = interval,
    data = summary_df
  )
  class(res) <- "summary.sb_beta"
  res
}

#' Print of summary of `sb_beta()` results
#' @rdname sb_beta_methods
#' @export
print.summary.sb_beta <- function(x, digits = 3, n = 10, ...) {
  cat("SelectBoost beta summary\n")
  if (!is.null(x$selector)) {
    cat("Selector: ", x$selector, "\n", sep = "")
  }
  if (!is.null(x$B)) {
    cat("Resamples per threshold: ", x$B, "\n", sep = "")
  }
  if (!is.null(x$interval)) {
    cat("Interval mode: ", x$interval, "\n", sep = "")
  }
  if (!is.null(x$c0)) {
    cat("c0 grid: ", paste(format(round(x$c0, digits), trim = TRUE), collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$steps.seq) && length(x$steps.seq)) {
    cat("Inner thresholds: ", paste(format(round(x$steps.seq, digits), trim = TRUE), collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$data) && nrow(x$data)) {
    display_n <- min(n, nrow(x$data))
    cat("Top rows:\n")
    print(utils::head(x$data, display_n), digits = digits, ...)
  } else {
    cat("(no summary data)\n")
  }
  invisible(x)
}


#' Autoplot for `sb_beta()` results
#' @rdname sb_beta_methods
#' @export
autoplot.sb_beta <- function(object, variables = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.sb_beta().", call. = FALSE)
  }
  freq_mat <- as.matrix(object)
  c0_levels <- attr(object, "c0.seq", exact = TRUE)
  if (is.null(c0_levels)) {
    c0_levels <- seq_len(nrow(freq_mat))
  }
  var_names <- colnames(freq_mat)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(ncol(freq_mat)))
  }
  plot_df <- data.frame(
    c0 = rep(c0_levels, each = ncol(freq_mat)),
    variable = rep(var_names, times = nrow(freq_mat)),
    frequency = as.numeric(freq_mat)
  )
  if (!is.null(variables)) {
    keep <- intersect(variables, var_names)
    if (!length(keep)) {
      stop("None of the requested variables were found in the sb_beta object.", call. = FALSE)
    }
    plot_df <- subset(plot_df, plot_df$variable %in% keep)
    plot_df$variable <- factor(plot_df$variable, levels = keep)
  } else {
    plot_df$variable <- factor(plot_df$variable, levels = var_names)
  }
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$c0, y = .data$frequency, colour = .data$variable)) +
    ggplot2::geom_line(...) +
    ggplot2::labs(
      x = expression(c[0]),
      y = "Selection frequency",
      colour = "Variable"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_minimal()
}
