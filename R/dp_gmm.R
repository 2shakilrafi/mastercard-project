# Differentially private univariate Gaussian mixture (noisy EM) and synthetic
# series generation.
#
# Privacy model
# -------------
# * Unit of protection: one row of a store's weekly-change table
#   (add/remove-one-record adjacency).
# * Mechanism: Gaussian mechanism applied to the EM sufficient statistics,
#   accounted with zero-concentrated DP (zCDP). Per-iteration costs add up;
#   the total rho is converted to (epsilon, delta) with the standard bound
#   epsilon = rho + 2 * sqrt(rho * log(1 / delta)).
# * Public inputs (never learned from the data): number of components K,
#   clipping bound, initial parameters, number of iterations, variance range.
# * Everything computed from released values (synthetic series, simulated
#   paths) is post-processing and costs no additional budget.
#
# Re-running a release on the same data with a new seed spends new budget.

# --- zCDP accounting ---------------------------------------------------------

zcdp_rho_from_eps <- function(epsilon, delta) {
  stopifnot(epsilon > 0, delta > 0, delta < 1)
  L <- log(1 / delta)
  (sqrt(L + epsilon) - sqrt(L))^2
}

zcdp_eps_from_rho <- function(rho, delta) {
  stopifnot(rho >= 0, delta > 0, delta < 1)
  rho + 2 * sqrt(rho * log(1 / delta))
}

# Standard deviation of Gaussian noise that gives rho-zCDP for an L2
# sensitivity of `sens`.
gaussian_sd <- function(sens, rho) {
  stopifnot(rho > 0)
  sens / sqrt(2 * rho)
}

# --- Sufficient statistics ---------------------------------------------------

# L2 bound on one record's contribution to the statistics vector below.
# Each record contributes (r, r * x, r * x^2) with sum(r) = 1 and |x| <= 1,
# so its L2 norm is at most sqrt(1 + 1 + 1).
dp_sensitivity <- function() sqrt(3)

clip_scale <- function(x, bound) {
  pmin(pmax(x, -bound), bound) / bound
}

# Responsibilities (n x K) under the current parameters, via log-sum-exp.
responsibilities <- function(xs, pi, mu, sigma) {
  K <- length(pi)
  lp <- vapply(seq_len(K), function(k) {
    log(pi[k]) + dnorm(xs, mu[k], sigma[k], log = TRUE)
  }, numeric(length(xs)))
  lp <- matrix(lp, nrow = length(xs), ncol = K)
  m <- apply(lp, 1, max)
  w <- exp(lp - m)
  w / rowSums(w)
}

# Per-record statistics summed over the data: counts, first and second
# moments for each component. Length 3K.
dp_stats <- function(xs, pi, mu, sigma) {
  r <- responsibilities(xs, pi, mu, sigma)
  c(colSums(r), colSums(r * xs), colSums(r * xs^2))
}

# The only place the data touch the released model: statistics plus
# calibrated Gaussian noise.
release_stats <- function(xs, par, noise_sd) {
  f <- dp_stats(xs, par$pi, par$mu, par$sigma)
  f + rnorm(length(f), 0, noise_sd)
}

# --- Fitting -----------------------------------------------------------------

default_init <- function(K) {
  list(
    pi = rep(1 / K, K),
    mu = if (K == 1) 0 else seq(-0.5, 0.5, length.out = K),
    sigma = rep(0.3, K)
  )
}

# Fit a K-component Gaussian mixture to weekly changes `x` (dollars) with a
# total zCDP budget `rho`. `bound` is a public clipping bound in dollars.
dp_gmm_fit <- function(x, K = 2, bound, rho, iters = 5,
                       sigma_range = c(0.02, 1), n_floor = 1, seed = NULL) {
  stopifnot(K >= 1, bound > 0, rho > 0, iters >= 1)
  if (!is.null(seed)) set.seed(seed)
  xs <- clip_scale(as.numeric(x[!is.na(x)]), bound)
  par <- default_init(K)
  rho_iter <- rho / iters
  noise_sd <- gaussian_sd(dp_sensitivity(), rho_iter)

  for (t in seq_len(iters)) {
    f <- release_stats(xs, par, noise_sd)
    # Post-processing of the noisy statistics from here on.
    n_k <- pmax(f[seq_len(K)], n_floor)
    s_k <- f[K + seq_len(K)]
    q_k <- f[2 * K + seq_len(K)]
    mu <- pmin(pmax(s_k / n_k, -1), 1)
    v <- pmin(pmax(q_k / n_k - mu^2, sigma_range[1]^2), sigma_range[2]^2)
    par <- list(pi = n_k / sum(n_k), mu = mu, sigma = sqrt(v))
  }

  structure(list(
    pi = par$pi, mu = par$mu * bound, sigma = par$sigma * bound,
    K = K, bound = bound, rho = rho, iters = iters,
    rho_per_iter = rho_iter, noise_sd = noise_sd,
    sensitivity = dp_sensitivity()
  ), class = "dp_gmm")
}

# Non-private reference fit with the same clipping, K, and initialization.
# Used only to measure the utility cost of privacy, never released.
gmm_fit_nonprivate <- function(x, K = 2, bound, iters = 200,
                               sigma_range = c(0.02, 1)) {
  xs <- clip_scale(as.numeric(x[!is.na(x)]), bound)
  par <- default_init(K)
  for (t in seq_len(iters)) {
    r <- responsibilities(xs, par$pi, par$mu, par$sigma)
    n_k <- pmax(colSums(r), 1e-8)
    mu <- colSums(r * xs) / n_k
    v <- pmin(pmax(colSums(r * xs^2) / n_k - mu^2, sigma_range[1]^2),
              sigma_range[2]^2)
    par <- list(pi = n_k / sum(n_k), mu = mu, sigma = sqrt(v))
  }
  structure(list(pi = par$pi, mu = par$mu * bound, sigma = par$sigma * bound,
                 K = K, bound = bound, rho = Inf, iters = iters),
            class = "dp_gmm")
}

# Public, data-free reference model (the initialization). A private model is
# only useful if it beats this.
gmm_uninformed <- function(K = 2, bound) {
  par <- default_init(K)
  structure(list(pi = par$pi, mu = par$mu * bound, sigma = par$sigma * bound,
                 K = K, bound = bound, rho = 0, iters = 0), class = "dp_gmm")
}

# --- Level -------------------------------------------------------------------

# Release the yearly mean sales level of a store with the Gaussian mechanism.
# `y_mean` is the mean of `n_weeks` weekly values, each assumed to lie in
# [0, s_max]; removing or adding one week moves the mean by at most
# s_max / (n_weeks - 1).
dp_level <- function(y_mean, s_max, rho, n_weeks = 52) {
  stopifnot(s_max > 0, rho > 0, n_weeks > 1)
  y <- pmin(pmax(y_mean, 0), s_max)
  sens <- s_max / (n_weeks - 1)
  pmin(pmax(y + rnorm(1, 0, gaussian_sd(sens, rho)), 0), s_max)
}

# --- Synthetic data (pure post-processing; takes no data) --------------------

mixture_sd <- function(model) {
  m <- sum(model$pi * model$mu)
  sqrt(sum(model$pi * (model$sigma^2 + model$mu^2)) - m^2)
}

mixture_mean <- function(model) sum(model$pi * model$mu)

sample_deltas <- function(model, n) {
  comp <- sample.int(model$K, n, replace = TRUE, prob = model$pi)
  rnorm(n, model$mu[comp], model$sigma[comp])
}

synthesize_series <- function(model, level, n_weeks = 260) {
  delta <- sample_deltas(model, n_weeks)
  data.frame(week = seq_len(n_weeks), delta = delta, sales = level + cumsum(delta))
}

# --- Budget ledger -----------------------------------------------------------

dp_ledger_row <- function(store, model, level_rho, level_bound, epsilon, delta,
                          seed, level_released) {
  rho_total <- model$rho + level_rho
  data.frame(
    store = store, K = model$K, clip_bound = model$bound, iters = model$iters,
    rho_model = model$rho, rho_level = level_rho, rho_total = rho_total,
    epsilon_target = epsilon, delta = delta,
    epsilon_spent = zcdp_eps_from_rho(rho_total, delta),
    level_bound = level_bound, level_released = level_released, seed = seed
  )
}
