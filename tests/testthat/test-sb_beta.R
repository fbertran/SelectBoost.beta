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
  
  expect_s3_class(print(fit), "sb_beta")
  expect_s3_class(print(out), "summary.sb_beta")
  
  skip_if_not_installed("ggplot2")
  plt <- autoplot.sb_beta(fit)
  expect_s3_class(plt, "ggplot")
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
  
  set.seed(42)
  expect_no_error(sb_resample_groups(X_norm, groups, B = 4, seed = 99))
  set.seed(42)
  expect_no_error(sb_resample_groups(X_norm, groups, B = 4, seed = 99, use.parallel = TRUE))
})
