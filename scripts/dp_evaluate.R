# Utility of the private fit versus the privacy budget.
# This script reads the PRIVATE inputs and runs the mechanism many times to
# study its accuracy. It is an analysis of the mechanism, not a release: only
# the summary tables it writes (dp-results/) are meant to be shared.
source("scripts/common.R")
set.seed(20261117)

reps    <- 200
delta   <- 1e-5
eps_set <- c(0.5, 1, 2, 4, 8)
K_set   <- c(1, 2)
B_set   <- c(5e5, 1e6)
stores  <- 1:5
share_em <- 0.9   # remaining 10% of rho goes to the level release

iters_for <- function(K) if (K == 1) 1 else 5

rows <- list()
for (i in stores) {
  emp <- load_store_deltas(i); emp <- emp[!is.na(emp)]
  for (K in K_set) for (B in B_set) {
    base <- list(
      nonprivate = gmm_fit_nonprivate(emp, K, B),
      uninformed = gmm_uninformed(K, B)
    )
    for (nm in names(base)) {
      m <- replicate(50, utility_metrics(base[[nm]], emp))
      rows[[length(rows) + 1]] <- data.frame(
        store = i, K = K, bound = B, epsilon = if (nm == "nonprivate") Inf else 0,
        method = nm, w1 = mean(m["w1", ]), w1_sd = sd(m["w1", ]),
        mean_err = mean(m["mean_err", ]), sd_ratio = mean(m["sd_ratio", ]))
    }
    for (eps in eps_set) {
      rho <- zcdp_rho_from_eps(eps, delta)
      m <- replicate(reps, {
        fit <- dp_gmm_fit(emp, K, B, share_em * rho, iters = iters_for(K))
        utility_metrics(fit, emp)
      })
      rows[[length(rows) + 1]] <- data.frame(
        store = i, K = K, bound = B, epsilon = eps, method = "dp",
        w1 = mean(m["w1", ]), w1_sd = sd(m["w1", ]),
        mean_err = mean(m["mean_err", ]), sd_ratio = median(m["sd_ratio", ]))
    }
  }
  message("store ", i, " done")
}
res <- do.call(rbind, rows)
dir.create("dp-results", showWarnings = FALSE)
write.csv(res, "dp-results/utility.csv", row.names = FALSE)

# --- Audit of the generator and of the inputs (non-private, diagnostic) -------
summ <- load_store_summary()
aud <- do.call(rbind, lapply(stores, function(i) {
  d <- load_store_deltas(i); lv <- cumsum(replace(d, is.na(d), 0))
  emp <- d[!is.na(d)]
  lab <- summ[summ$store == paste("Store", i), ]
  match_row <- summ$store[which.min(abs(summ$level_sd - sd(lv)))]
  m <- gmm_fit_nonprivate(emp, 2, 1e6)
  sim_sd <- mean(replicate(300, sd(synthesize_series(m, 0, 52)$sales)))
  data.frame(
    store = i,
    summary_level_mean = lab$level_mean, summary_level_sd = lab$level_sd,
    recon_level_sd = sd(lv),
    sd_matches_labeled_row = abs(sd(lv) / lab$level_sd - 1) < 0.01,
    sd_closest_row = match_row,
    lag1_autocorr_real = cor(emp[-length(emp)], emp[-1]),
    sim_52wk_sd_nonprivate_K2 = sim_sd,
    old_shift_constant = 2425553)
}))
write.csv(aud, "dp-results/audit.csv", row.names = FALSE)
print(aud)
