# Shared helpers for the DP scripts. Run from the repository root.
source("R/dp_gmm.R")

load_store_deltas <- function(i, dir = "raw-csv") {
  d <- read.csv(file.path(dir, sprintf("store-%d.csv", i)))
  as.numeric(d$difference_from_previous_week)
}

load_store_summary <- function(dir = "raw-csv") {
  s <- read.csv(file.path(dir, "all-store-avg-+-stddev.csv"), check.names = FALSE)
  names(s) <- c("store", "level_mean", "level_sd")
  s
}

# Wasserstein-1 distance between a model and an empirical sample, computed on
# a quantile grid, in units of the empirical standard deviation.
w1_to_empirical <- function(model, emp, m = 2000) {
  p <- (seq_len(m) - 0.5) / m
  q_emp <- quantile(emp, p, type = 7, names = FALSE)
  q_mod <- sort(sample_deltas(model, m))
  mean(abs(q_mod - q_emp)) / sd(emp)
}

utility_metrics <- function(model, emp) {
  s <- sd(emp)
  c(w1 = w1_to_empirical(model, emp),
    mean_err = abs(mixture_mean(model) - mean(emp)) / s,
    sd_ratio = mixture_sd(model) / s)
}
