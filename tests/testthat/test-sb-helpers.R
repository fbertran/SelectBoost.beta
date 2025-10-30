test_that("sb_normalize centres and scales columns", {
  set.seed(123)
  X <- matrix(rnorm(200), ncol = 4)
  X_norm <- sb_normalize(X)
  expect_equal(colMeans(X_norm), rep(0, ncol(X)), tolerance = 1e-10)
  expect_equal(sqrt(colSums(X_norm^2)), rep(1, ncol(X)), tolerance = 1e-10)
  expect_null(colnames(X_norm))
  expect_equal(attr(X_norm, "center"), colMeans(X))
  expect_equal(attr(X_norm, "scale"), sqrt(colSums(sweep(X,2,colMeans(X),`-`)^2)))
})

test_that("sb_normalize handles custom centering/scaling", {
  X <- matrix(rep(1:5, each = 3), ncol = 3)
  cn <- c(2, 3, 4)
  sc <- c(1, 2, 3)
  X_norm <- sb_normalize(X, center = cn, scale = sc)
  expect_equal(colMeans(X_norm), colMeans(scale(X,center=cn,scale=sc)), tolerance = 1e-10)
  expect_equal(attr(X_norm, "center"), rep(cn, length.out = ncol(X)))
  expect_equal(attr(X_norm, "scale"), rep(sc, length.out = ncol(X)))
})

test_that("sb_resample_groups generates correlated draws", {
  set.seed(99)
  X <- matrix(rnorm(1000), ncol = 4)
  X[, 2] <- X[, 1] * 0.8 + rnorm(nrow(X), sd = 0.1)
  groups <- c(1, 1, 2, 2)
  mat <- sb_resample_groups(X, groups, B = 2000, seed = 42)
  expect_equal(length(mat), 2000)
  expect_equal(dim(mat[[1]]),c(250, 4))
  expect_null(colnames(mat[[1]]))
  expect_gt(cor(mat[[1]][, 1], mat[[1]][, 2]), 0.4)
  expect_lt(abs(cor(mat[[1]][, 1], mat[[1]][, 3])), 0.2)
  expect_equal(attr(mat[[1]], "groups"), groups)
})

test_that("sb_normalize centres and scales columns", {
  set.seed(123)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  X_norm <- sb_normalize(X)
  zero = rep(0, ncol(X))
  names(zero) <- paste0("x", seq_len(ncol(X)))
  expect_equal(colMeans(X_norm), zero, tolerance = 1e-8)
  ones = rep(1, ncol(X))
  names(ones) <- paste0("x", seq_len(ncol(X)))
  expect_equal(sqrt(colSums(X_norm^2)), ones, tolerance = 1e-8)
})

test_that("sb_group_variables builds correlation groups", {
  corr_mat <- matrix(
    c(1, 0.8, 0.65, 0.2,
      0.8, 1, 0.61, 0.15,
      0.65, 0.61, 1, 0.05,
      0.2, 0.15, 0.05, 1),
    nrow = 4, byrow = TRUE
  )
  groups <- sb_group_variables(corr_mat, c0 = 0.6)
  
  expect_type(groups, "list")
  expect_true(all(vapply(groups, is.integer, logical(1))))
  expect_equal(groups[[1]], c(1L, 2L, 3L))
  expect_equal(groups[[4]], 4L)
})

test_that("sb_resample_groups generates correlated draws", {
  set.seed(456)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  X_norm <- sb_normalize(X)
  groups <- list(1:2, 1:2, 3L, 4L)
  
  draws <- sb_resample_groups(X_norm, groups, B = 3, seed = 789)
  expect_length(draws, 3)
  lapply(draws, function(mat) {
    expect_equal(dim(mat), dim(X_norm))
    expect_equal(rownames(mat), NULL)
    expect_equal(colnames(mat), colnames(X_norm))
    zeros <- c(0, 0)
    names(zeros) <- paste0("x", seq_len(length(zeros)))
    expect_equal(colMeans(mat[, 1:2, drop = FALSE]), zeros, tolerance = 1e-7)
    ones <- c(1, 1)
    names(ones) <- paste0("x", seq_len(length(ones)))
    expect_equal(sqrt(colSums(mat[, 1:2, drop = FALSE]^2)), ones, tolerance = 1e-7)
  })
})

test_that("sb_selection_frequency handles glmnet and lars conventions", {
  mat <- matrix(c(
    0.2, 0.0,
    0.9, 0.1,
    -0.05, 0.2
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- c("(Intercept)", "x1", "x2")
  colnames(mat) <- c("sim1", "sim2")
  
  freq_glmnet <- sb_selection_frequency(mat, version = "glmnet", threshold = 0.15)
  expect_equal(freq_glmnet, c(x1 = 0.5, x2 = 0.5))
  
  freq_lars <- sb_selection_frequency(mat[-1, , drop = FALSE], version = "lars", threshold = 0.15)
  expect_equal(freq_lars, c(x1 = 0.5, x2 = 0.5))
})

test_that("sb_beta accepts custom selectors", {
  set.seed(321)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  Y <- stats::runif(10, min = 0.2, max = 0.8)
  
  toy_selector <- function(X, Y, ...) {
    nm <- colnames(X)
    res <- rep(0, length(nm) + 1)
    names(res) <- c("(Intercept)", nm)
    res[2] <- 1
    res
  }
  
  res <- sb_beta(
    X, Y,
    selector = toy_selector,
    B = 4,
    step.num = 0.5,
    steps.seq = c(0.75, 0.5),
    version = "glmnet",
    threshold = 0.5
  )
  
  expect_s3_class(res, "sb_beta")
  expect_equal(attr(res, "B"), 4)
  expect_equal(attr(res, "steps.seq"), c(0.75, 0.5))
  expect_equal(attr(res, "selector"), "selector")
  expect_equal(rownames(res), sprintf("c0 = %.3f", c(1, 0.75, 0.5, 0)))
  expect_true(all(res[, "x1"] == 1))
  expect_true(all(res[, colnames(X)[-1]] == 0))
  
  attr(toy_selector,"fun.name") <- "toy_selector"
  res <- sb_beta(
    X, Y,
    selector = toy_selector,
    B = 4,
    step.num = 0.5,
    steps.seq = c(0.75, 0.5),
    version = "glmnet",
    threshold = 0.5
  )
  
  expect_s3_class(res, "sb_beta")
  expect_equal(attr(res, "B"), 4)
  expect_equal(attr(res, "steps.seq"), c(0.75, 0.5))
  expect_equal(attr(res, "selector"), "toy_selector")
  expect_equal(rownames(res), sprintf("c0 = %.3f", c(1, 0.75, 0.5, 0)))
  expect_true(all(res[, "x1"] == 1))
  expect_true(all(res[, colnames(X)[-1]] == 0))
})

