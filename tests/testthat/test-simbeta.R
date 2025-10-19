test_that("simulation_DATA.beta produces valid intervals and missingness", {
  set.seed(123)
  sim <- simulation_DATA.beta(
    n = 120, p = 8, s = 3, beta_size = 1.0,
    corr = "ar1", rho = 0.2,
    mechanism = "mixed", mix_prob = 0.5,
    delta = 0.04, alpha = 0.1, na_rate = 0.15, na_side = "random"
  )
  expect_equal(dim(sim$X), c(120,8))
  expect_true(all(sim$Y > 0 & sim$Y < 1))
  ok <- is.finite(sim$Y_low) & is.finite(sim$Y_high)
  expect_true(mean(!ok) >= 0.1)
  expect_true(all(sim$Y_low[ok] >= 0 & sim$Y_low[ok] <= 1))
  expect_true(all(sim$Y_high[ok] >= 0 & sim$Y_high[ok] <= 1))
  expect_true(all(sim$Y_low[ok] <= sim$Y_high[ok]))
})
