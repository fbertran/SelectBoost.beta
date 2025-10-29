test_that("betareg stepwise selectors support phi-submodel selection", {
  skip_if_not_installed("betareg")
  
  set.seed(11)
  n <- 80
  X_mean <- matrix(rnorm(n * 2), n, 2)
  colnames(X_mean) <- c("x1", "x2")
  X_phi <- matrix(rnorm(n * 2), n, 2)
  colnames(X_phi) <- c("z1", "z2")
  mu_lin <- 0.5 + 1.2 * X_mean[, 1] - 0.8 * X_mean[, 2]
  mu <- plogis(mu_lin)
  phi_lin <- 1 + 1.5 * X_phi[, 1]
  phi <- exp(phi_lin)
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  
  coefs <- betareg_step_aic(
    X_mean,
    y,
    X_phi = X_phi,
    direction_phi = "forward",
    max_steps = 6
  )
  
  phi_terms <- coefs[grepl("^phi\\|", names(coefs))]
  expect_true(any(names(phi_terms) == "phi|z1" && abs(phi_terms["phi|z1"]) > 0))
  expect_true("(Intercept)" %in% names(coefs))
})

test_that("betareg stepwise selectors honour observation weights", {
  skip_if_not_installed("betareg")
  
  set.seed(12)
  n <- 60
  X <- matrix(rnorm(n * 2), n, 2)
  colnames(X) <- c("x1", "x2")
  mu <- plogis(0.5 + 0.6 * X[, 1])
  y <- rbeta(n, mu * 15, (1 - mu) * 15)
  w <- runif(n, 0.5, 2)
  
  weighted_fit <- betareg_step_aic(X, y, weights = w, max_steps = 0)
  unweighted_fit <- betareg_step_aic(X, y, max_steps = 0)
  
  ref <- betareg::betareg(
    y ~ 1,
    data = data.frame(y = y, X),
    weights = w
  )
  ref_intercept <- unname(coef(ref)$mean["(Intercept)"])
  
  expect_equal(weighted_fit["(Intercept)"], ref_intercept, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(weighted_fit["(Intercept)"], unweighted_fit["(Intercept)"])))
})
