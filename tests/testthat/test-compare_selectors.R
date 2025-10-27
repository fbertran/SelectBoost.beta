test_that("compare_selectors_single drops incomplete cases", {
  set.seed(2024)
  X <- matrix(rnorm(300), nrow = 60, ncol = 5)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  Y <- plogis(X[, 1] - 0.3 * X[, 2])
  Y <- rbeta(nrow(X), Y * 15, (1 - Y) * 15)
  Y[c(1, 5)] <- NA_real_
  X[c(2, 3), 4] <- NA_real_
  
  expect_error(
    compare_selectors_single(X, Y, include_enet = FALSE),
    NA
  )
  
  expect_error(
    compare_selectors_bootstrap(X, Y, B = 4, include_enet = FALSE, seed = 99),
    NA
  )
})

test_that("compare selectors guard against fully missing data", {
  X <- matrix(rnorm(40), nrow = 20, ncol = 2)
  Y <- rep(NA_real_, nrow(X))
  
  expect_error(
    compare_selectors_single(X, Y, include_enet = FALSE),
    "No complete cases"
  )
  
  expect_error(
    compare_selectors_bootstrap(X, Y, include_enet = FALSE),
    "No complete cases"
  )
})
