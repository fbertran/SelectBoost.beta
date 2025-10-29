test_that("sb_beta performs correlated resampling and tallies frequencies", {
  skip_if_not_installed("betareg")
  
  set.seed(2025)
  sim <- simulation_DATA.beta(n = 40, p = 4, s = 2, beta_size = 0.7)
  fit <- sb_beta(sim$X, sim$Y, B = 5, step.num = 0.5, version = "glmnet", threshold = 1e-6)
  
  expect_true("sb_beta" %in% class(fit))
  expect_equal(ncol(fit), ncol(sim$X))
  expect_equal(attr(fit, "B"), 5)
  expect_identical(attr(fit, "selector"), "betareg_step_aic")
  expect_true(all(is.finite(fit)))
  expect_true(all(fit >= 0 & fit <= 1))
  expect_true(any(fit[nrow(fit), ] > 0))
  diag_list <- attr(fit, "resample_diagnostics")
  expect_true(is.list(diag_list))
  expect_equal(length(diag_list), nrow(fit))
  if (length(diag_list)) {
    expect_s3_class(diag_list[[length(diag_list)]], "data.frame")
  }
})

test_that("sb_beta summary and plotting helpers work", {
  skip_if_not_installed("betareg")
  
  set.seed(123)
  sim <- simulation_DATA.beta(n = 30, p = 3, s = 2, beta_size = 0.4)
  fit <- sb_beta(sim$X, sim$Y, B = 3, step.num = 0.5, version = "glmnet")
  
  out <- summary(fit)
  expect_s3_class(out, "summary.sb_beta")
  expect_equal(out$selector, "betareg_step_aic")
  expect_equal(out$B, 3)
  expect_true(is.data.frame(out$data))
  expect_true(all(out$data$frequency >= 0 & out$data$frequency <= 1))
  expect_equal(attr(fit, "interval"), "none")
  
  expect_s3_class(print(fit), "sb_beta")
  expect_s3_class(print(out), "summary.sb_beta")
  
  skip_if_not_installed("ggplot2")
  plt <- autoplot.sb_beta(fit)
  expect_s3_class(plt, "ggplot")
})

test_that("sb_resample_groups exposes diagnostics and caching", {
  set.seed(2027)
  X <- matrix(rnorm(80), nrow = 20)
  X_norm <- sb_normalize(X)
  corr <- sb_compute_corr(X_norm)
  groups <- sb_group_variables(corr, c0 = 0.2)
  
  cache_env <- new.env(parent = emptyenv())
  first <- sb_resample_groups(X_norm, groups, B = 3, cache = cache_env)
  second <- sb_resample_groups(X_norm, groups, B = 3, cache = cache_env)
  
  expect_equal(length(first), 3L)
  expect_equal(length(second), 3L)
  for (b in seq_along(first)) {
    expect_equal(first[[b]], second[[b]])
  }
  diag_first <- attr(first, "diagnostics")
  diag_second <- attr(second, "diagnostics")
  expect_true(is.data.frame(diag_first))
  expect_true(is.data.frame(diag_second))
  expect_true(any(diag_first$regenerated > 0 | diag_first$size < 2))
  if (nrow(diag_second)) {
    expect_true(all(diag_second$regenerated[diag_second$size >= 2] == 0))
    expect_true(all(diag_second$cached[diag_second$size >= 2]))
  }
})

test_that("sb_resample_groups can reuse future.apply for parallelism", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  skip_on_cran()
  
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential()), add = TRUE)
  
  set.seed(2026)
  X <- matrix(rnorm(60), nrow = 20)
  X_norm <- sb_normalize(X)
  corr <- sb_compute_corr(X_norm)
  groups <- sb_group_variables(corr, c0 = 0.2)
  
  cache_env <- new.env(parent = emptyenv())
  set.seed(42)
  seq_res <- sb_resample_groups(X_norm, groups, B = 4, cache = cache_env)
  set.seed(42)
  expect_no_error(sb_resample_groups(X_norm, groups, B = 4, seed = 99, use.parallel = TRUE))
  
  par_res <- sb_resample_groups(X_norm, groups, B = 4, use.parallel = TRUE, cache = cache_env)
  
  for (b in seq_along(seq_res)) {
    expect_equal(par_res[[b]], seq_res[[b]])
  }
  expect_equal(attr(seq_res, "diagnostics"), attr(par_res, "diagnostics"))
})

test_that("sb_beta handles interval outcomes", {
  skip_if_not_installed("betareg")
  
  set.seed(101)
  sim <- simulation_DATA.beta(n = 30, p = 3, s = 1, beta_size = 0.4)
  y_low <- pmax(sim$Y - 0.05, 0)
  y_high <- pmin(sim$Y + 0.05, 1)
  fit <- sb_beta(sim$X, Y_low = y_low, Y_high = y_high, interval = "uniform", B = 4, step.num = 0.5)
  expect_equal(attr(fit, "interval"), "uniform")
  diag_list <- attr(fit, "resample_diagnostics")
  expect_equal(length(diag_list), nrow(fit))
  expect_true(all(vapply(diag_list, is.data.frame, logical(1))))
})
