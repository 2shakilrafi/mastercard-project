library(testthat)
source(file.path("..", "..", "R", "dp_gmm.R"))

# small helper so the tests need no extra package
gtools_rdirichlet <- function(K) {
  g <- rgamma(K, 1)
  g / sum(g)
}

test_that("zCDP <-> (epsilon, delta) conversion round-trips", {
  for (eps in c(0.1, 1, 4, 10)) {
    for (delta in c(1e-3, 1e-5, 1e-8)) {
      rho <- zcdp_rho_from_eps(eps, delta)
      expect_equal(zcdp_eps_from_rho(rho, delta), eps, tolerance = 1e-9)
    }
  }
})

test_that("gaussian_sd gives exactly the requested rho", {
  sens <- 1.7
  rho <- 0.05
  sd <- gaussian_sd(sens, rho)
  expect_equal(sens^2 / (2 * sd^2), rho)
})

test_that("per-record statistics respect the stated L2 sensitivity", {
  set.seed(1)
  worst <- 0
  for (trial in 1:3000) {
    K <- sample(1:4, 1)
    pi <- as.numeric(gtools_rdirichlet(K))
    mu <- runif(K, -1, 1)
    sigma <- runif(K, 0.02, 1)
    n <- sample(0:30, 1)
    base <- if (n > 0) runif(n, -1, 1) else numeric(0)
    # include boundary values, where the bound is tightest
    extra <- sample(c(-1, 1, runif(1, -1, 1)), 1)
    d <- dp_stats(c(base, extra), pi, mu, sigma) -
      (if (n > 0) dp_stats(base, pi, mu, sigma) else rep(0, 3 * K))
    worst <- max(worst, sqrt(sum(d^2)))
  }
  expect_lte(worst, dp_sensitivity() + 1e-9)
})

test_that("noise added to the statistics has the calibrated standard deviation", {
  x <- c(-2e5, 1e5, 3e5, -4e4, 2.5e5, 0)
  bound <- 1e6
  rho <- 0.2
  xs <- clip_scale(x, bound)
  par <- default_init(2)
  sd_theory <- gaussian_sd(dp_sensitivity(), rho)
  set.seed(2)
  clean <- dp_stats(xs, par$pi, par$mu, par$sigma)
  noise <- replicate(20000, release_stats(xs, par, sd_theory) - clean)
  expect_equal(apply(noise, 1, sd), rep(sd_theory, 6), tolerance = 0.03)
  expect_equal(abs(rowMeans(noise)), rep(0, 6), tolerance = 0.1 * sd_theory)
  # the fit uses exactly this calibration, split evenly over iterations
  fit <- dp_gmm_fit(x, K = 2, bound = bound, rho = rho, iters = 4, seed = 3)
  expect_equal(fit$noise_sd, gaussian_sd(dp_sensitivity(), rho / 4))
  expect_equal(fit$rho_per_iter * fit$iters, rho)
})

test_that("with a huge budget the fit recovers a well-separated mixture", {
  set.seed(4)
  bound <- 1
  n <- 6000
  comp <- rbinom(n, 1, 0.3) + 1
  x <- rnorm(n, c(-0.4, 0.4)[comp], c(0.1, 0.15)[comp])
  fit <- dp_gmm_fit(x, K = 2, bound = bound, rho = 1e12, iters = 60, seed = 5)
  o <- order(fit$mu)
  expect_equal(fit$mu[o], c(-0.4, 0.4), tolerance = 0.03)
  expect_equal(fit$sigma[o], c(0.1, 0.15), tolerance = 0.03)
  expect_equal(fit$pi[o], c(0.7, 0.3), tolerance = 0.03)
})

test_that("fit output is valid and reproducible under a seed", {
  x <- rnorm(51, 0, 3e5)
  a <- dp_gmm_fit(x, K = 2, bound = 1e6, rho = 0.1, seed = 10)
  b <- dp_gmm_fit(x, K = 2, bound = 1e6, rho = 0.1, seed = 10)
  expect_equal(a$mu, b$mu)
  expect_equal(sum(a$pi), 1)
  expect_true(all(a$sigma >= 0.02 * 1e6 - 1e-6 & a$sigma <= 1e6 + 1e-6))
  expect_true(all(abs(a$mu) <= 1e6 + 1e-6))
})

test_that("values beyond the public bound are clipped, not trusted", {
  x1 <- c(1e5, -2e5, 3e5)
  x2 <- c(1e5, -2e5, 9e9)  # outlier far outside the bound
  x3 <- c(1e5, -2e5, 1e6)  # exactly at the bound
  expect_equal(clip_scale(x2, 1e6), clip_scale(x3, 1e6))
  expect_equal(dp_stats(clip_scale(x2, 1e6), 1, 0, 0.3),
               dp_stats(clip_scale(x3, 1e6), 1, 0, 0.3))
})

test_that("level release is clipped to its public range and calibrated", {
  set.seed(6)
  v <- replicate(5000, dp_level(2.4e6, 5e6, rho = 0.05))
  expect_true(all(v >= 0 & v <= 5e6))
  sd_theory <- gaussian_sd(5e6 / 51, 0.05)
  expect_equal(sd(v), sd_theory, tolerance = 0.05)
})

test_that("synthetic generation needs no data and returns the right shape", {
  m <- gmm_uninformed(2, 1e6)
  s <- synthesize_series(m, level = 2e6, n_weeks = 40)
  expect_equal(nrow(s), 40)
  expect_equal(s$sales, 2e6 + cumsum(s$delta))
})

test_that("ledger totals match the budget that was requested", {
  eps <- 4
  delta <- 1e-5
  rho <- zcdp_rho_from_eps(eps, delta)
  model <- dp_gmm_fit(rnorm(30, 0, 2e5), K = 2, bound = 1e6, rho = 0.9 * rho,
                      seed = 1)
  row <- dp_ledger_row("Store 1", model, 0.1 * rho, 5e6, eps, delta, 1, 2.4e6)
  expect_equal(row$rho_total, rho)
  expect_equal(row$epsilon_spent, eps, tolerance = 1e-9)
})
