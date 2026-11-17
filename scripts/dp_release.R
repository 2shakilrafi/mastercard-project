# One differentially private release per store. Reads PRIVATE inputs
# (raw-csv/) and writes only privatized outputs to dp-synth-data/.
# Re-running with a different seed spends additional privacy budget on the
# same data; keep the ledger and do not release more than one run per store.
source("scripts/common.R")

config <- list(
  epsilon = 8, delta = 1e-5,
  K = 1,                 # fixed in advance of the fit; see dp-analysis for why
  clip_bound = 5e5,      # public bound on |weekly change|, dollars
  level_bound = 5e6,     # public upper bound on yearly mean level, dollars
  iters = 1, n_weeks_out = 52, share_model = 0.9
)
stopifnot(config$iters == (if (config$K == 1) 1 else config$iters))

rho_total <- zcdp_rho_from_eps(config$epsilon, config$delta)
rho_model <- config$share_model * rho_total
rho_level <- rho_total - rho_model
summ <- load_store_summary()
sd_check <- read.csv("dp-results/audit.csv")

dir.create("dp-synth-data", showWarnings = FALSE)
ledger <- list(); params <- list()
for (i in 1:5) {
  seed <- 20261117 + i
  set.seed(seed)
  d <- load_store_deltas(i)
  model <- dp_gmm_fit(d, K = config$K, bound = config$clip_bound,
                      rho = rho_model, iters = config$iters)
  true_mean <- summ$level_mean[summ$store == paste("Store", i)]
  level <- dp_level(true_mean, config$level_bound, rho_level)
  series <- synthesize_series(model, level, config$n_weeks_out)
  write.csv(series, sprintf("dp-synth-data/store-%d-dp-synth-sales.csv", i),
            row.names = FALSE)
  row <- dp_ledger_row(sprintf("Store %d", i), model, rho_level,
                       config$level_bound, config$epsilon, config$delta, seed,
                       level)
  row$level_source_check <- if (sd_check$sd_matches_labeled_row[i]) "ok" else
    paste("sd matches", sd_check$sd_closest_row[i], "not the labeled row")
  ledger[[i]] <- row
  params[[i]] <- data.frame(store = sprintf("Store %d", i), component = seq_len(model$K),
                            weight = model$pi, mean = model$mu, sd = model$sigma)
}
ledger <- do.call(rbind, ledger)
ledger$level_released <- round(ledger$level_released)
write.csv(ledger, "dp-synth-data/privacy-ledger.csv", row.names = FALSE)
write.csv(do.call(rbind, params), "dp-synth-data/dp-model-params.csv", row.names = FALSE)
print(ledger[, c("store","K","clip_bound","rho_total","epsilon_spent","delta","level_released","level_source_check")])
print(do.call(rbind, params))
