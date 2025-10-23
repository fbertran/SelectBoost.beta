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
