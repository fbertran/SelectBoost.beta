#' Simulate **interval** Beta-regression data (flexible)
# #' @inheritParams stats::rnorm
#' @param n,p Sample size and number of predictors.
#' @param s Number of active (nonzero) coefficients.
#' @param beta_size Scalar (alternating ±) or numeric vector of length ≥ s.
#' @param a0 Intercept (logit scale).
#' @param X_dist Distribution for X: `"gaussian"`, `"t"`, or `"bernoulli"`.
#' @param corr Correlation structure: `"indep"`, `"ar1"`, or `"block"`.
#' @param rho AR(1) correlation or within-block correlation.
#' @param block_size Block size when `corr = "block"` (default 5).
#' @param df Degrees of freedom for `X_dist = "t"` (default 5).
#' @param prob Success prob for `X_dist = "bernoulli"` (default 0.5).
#' @param active_idx Optional integer vector of active feature indices (length s). If NULL, uses 1:s.
#' @param phi Precision parameter: scalar, length-n vector, or function `(mu, X) -> length-n`.
#' @param mechanism Interval mechanism per row: `"jitter"`, `"quantile"`, or `"mixed"`.
#' @param mix_prob Probability of jitter when `mechanism = "mixed"`.
#' @param delta Symmetric jitter half-width (scalar / vector / function).
#' @param delta_low,delta_high Asymmetric jitter widths (override `delta` if set).
#' @param alpha Miscoverage for quantile intervals (scalar / vector / function).
#' @param alpha_low,alpha_high Asymmetric miscoverage (override `alpha` if set).
#' @param na_rate Fraction of rows with a missing bound (default 0).
#' @param na_side Which bound to drop: `"left"`, `"right"`, or `"random"`.
#' @param centerX,scaleX Whether to center/scale X before returning.
#' @param seed RNG seed.
#' @return list with `X, Y, Y_low, Y_high, mu, beta, a0, phi, info, settings`.
#' @export
simulation_DATA.beta <- function(
  n, p, s = min(5L, p), beta_size = 1, a0 = 0,
  X_dist = c("gaussian","t","bernoulli"),
  corr = c("indep","ar1","block"),
  rho = 0.0, block_size = 5L, df = 5, prob = 0.5,
  active_idx = NULL,
  phi = 20,
  mechanism = c("jitter","quantile","mixed"),
  mix_prob = 0.5,
  delta = 0.05,
  delta_low = NULL, delta_high = NULL,
  alpha = 0.1,
  alpha_low = NULL, alpha_high = NULL,
  na_rate = 0.0, na_side = c("random","left","right"),
  centerX = FALSE, scaleX = FALSE,
  seed = NULL
) {
  X_dist <- match.arg(X_dist)
  corr   <- match.arg(corr)
  mechanism <- match.arg(mechanism)
  na_side <- match.arg(na_side)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n > 0, p > 0, s >= 0, s <= p)

  # --- Design matrix X ---
  make_Sigma <- function(p, rho, block_size) {
    if (corr == "ar1") {
      idx <- seq_len(p)
      outer(idx, idx, function(i,j) rho^abs(i-j))
    } else if (corr == "block") {
      k <- ceiling(p / block_size)
      S <- diag(p)
      for (b in seq_len(k)) {
        a <- (b-1)*block_size + 1
        z <- min(b*block_size, p)
        S[a:z, a:z] <- rho
        diag(S)[a:z] <- 1
      }
      S
    } else diag(p)
  }
  Sigma <- make_Sigma(p, rho, block_size)

  if (corr == "indep") {
    Z <- switch(X_dist,
      gaussian = matrix(rnorm(n*p), n, p),
      t        = matrix(rt(n*p, df = df), n, p),
      bernoulli= matrix(rbinom(n*p, 1, prob), n, p)
    )
  } else {
    Z0 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    if (X_dist == "gaussian") {
      Z <- Z0
    } else if (X_dist == "t") {
      g <- rgamma(n, shape = df/2, rate = df/2)
      Z <- Z0 / sqrt(g)
    } else {
      thr <- qnorm(prob)
      Z <- (Z0 > thr) * 1
    }
  }
  colnames(Z) <- paste0("x", seq_len(p))
  X <- Z

  if (centerX) X <- scale(X, center = TRUE, scale = FALSE)
  if (scaleX)  X <- scale(X, center = FALSE, scale = TRUE)

  # --- Coefficients ---
  if (is.null(active_idx)) active_idx <- seq_len(s)
  stopifnot(length(active_idx) == s, max(active_idx) <= p)
  beta <- numeric(p)
  if (length(beta_size) == 1L) {
    beta[active_idx] <- rep(c(1,-1), length.out = s) * abs(beta_size)
  } else {
    if (length(beta_size) < s) stop("beta_size length < s")
    beta[active_idx] <- beta_size[seq_len(s)]
  }

  # --- Latent mean and precision ---
  eta <- as.vector(a0 + X %*% beta)
  mu  <- plogis(eta)

  if (is.function(phi)) {
    phi_vec <- as.numeric(phi(mu, X))
    if (length(phi_vec) != n) stop("phi(mu, X) must return length-n vector")
  } else if (length(phi) == 1L) {
    phi_vec <- rep(as.numeric(phi), n)
  } else {
    phi_vec <- as.numeric(phi); stopifnot(length(phi_vec) == n)
  }

  a <- mu * phi_vec
  b <- (1 - mu) * phi_vec
  Y <- rbeta(n, a, b)

  # --- helper to evaluate possibly functional args ---
  as_vec <- function(obj, default) {
    if (is.null(obj)) return(rep(default, n))
    if (is.function(obj)) {
      out <- as.numeric(obj(mu, X)); if (length(out) != n) stop("function must return length-n vector")
      return(out)
    }
    if (length(obj) == 1L) return(rep(as.numeric(obj), n))
    out <- as.numeric(obj); if (length(out) != n) stop("must be scalar or length-n")
    out
  }

  # --- per-row mechanism ---
  if (mechanism == "mixed") {
    mp <- as_vec(mix_prob, 0.5)
    mech <- ifelse(runif(n) < mp, "jitter", "quantile")
  } else mech <- rep(mechanism, n)

  # --- jitter widths ---
  base_delta <- as_vec(delta, 0.05)
  dL <- if (!is.null(delta_low)) as_vec(delta_low, 0.05) else base_delta
  dH <- if (!is.null(delta_high)) as_vec(delta_high, 0.05) else base_delta

  # --- quantile levels (asymmetric allowed) ---
  a_base <- as_vec(alpha, 0.1)
  aL <- if (!is.null(alpha_low)) as_vec(alpha_low, 0.1) else a_base
  aH <- if (!is.null(alpha_high)) as_vec(alpha_high, 0.1) else a_base
  aL <- pmin(pmax(aL, 1e-6), 0.99); aH <- pmin(pmax(aH, 1e-6), 0.99)

  # --- build intervals ---
  Y_low  <- numeric(n); Y_high <- numeric(n)
  for (i in seq_len(n)) {
    if (mech[i] == "jitter") {
      lo <- Y[i] - dL[i]; hi <- Y[i] + dH[i]
      Y_low[i]  <- max(0, lo)
      Y_high[i] <- min(1, hi)
    } else {
      al <- aL[i] / 2; ah <- aH[i] / 2
      lo <- qbeta(al, a[i], b[i]); hi <- qbeta(1 - ah, a[i], b[i])
      Y_low[i]  <- max(0, lo)
      Y_high[i] <- min(1, hi)
    }
  }
  swap <- Y_low > Y_high
  if (any(swap)) { tmp <- Y_low[swap]; Y_low[swap] <- Y_high[swap]; Y_high[swap] <- tmp }

  # --- optional missing bounds ---
  if (na_rate > 0) {
    m <- ceiling(na_rate * n)
    rows <- sample.int(n, m)
    if (na_side == "left") {
      Y_low[rows] <- NA_real_
    } else if (na_side == "right") {
      Y_high[rows] <- NA_real_
    } else {
      pick_left <- runif(m) < 0.5
      Y_low[rows[pick_left]]  <- NA_real_
      Y_high[rows[!pick_left]] <- NA_real_
    }
  }

  info <- list(mechanism = mech, delta_low = dL, delta_high = dH,
               alpha_low = aL, alpha_high = aH)

  settings <- list(
    n=n,p=p,s=s,beta_size=beta_size,a0=a0,X_dist=X_dist,corr=corr,rho=rho,
    block_size=block_size,df=df,prob=prob,active_idx=active_idx,
    mechanism=mechanism,mix_prob=mix_prob,centerX=centerX,scaleX=scaleX,seed=seed
  )

  list(X = X, Y = Y, Y_low = Y_low, Y_high = Y_high,
       mu = mu, beta = beta, a0 = a0, phi = phi_vec,
       info = info, settings = settings)
}
