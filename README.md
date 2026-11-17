# Gaussian Mixture Modeling Project

An exploratory University of Arkansas and Mastercard project on generating synthetic store-level sales series from weekly summary statistics, using Gaussian mixture models (GMMs) in R.

There are two pipelines here, and only one of them is privacy preserving:

| Pipeline | Files | Privacy |
| --- | --- | --- |
| **Differentially private** (current) | `R/dp_gmm.R`, `scripts/`, `dp-analysis.qmd`, `dp-synth-data/`, `tests/` | $(\varepsilon, \delta)$-DP, $\varepsilon = 8$, $\delta = 10^{-5}$ per store, one release per store |
| **Non-private baseline** (original) | `analysis.qmd`, `analysis.html`, `recon/`, `synth-data/` | None. Do not release its outputs. |

Read [`dp-analysis.qmd`](dp-analysis.qmd) (rendered as `dp-analysis.html`) for the method, the privacy budget, the accuracy cost, and the problems found in the original analysis.

## What the private pipeline guarantees

- The unit of protection is one row of a store's weekly-change table.
- The model is fit with noisy EM: Gaussian noise is added to the sufficient statistics, calibrated to a proven sensitivity bound and tracked with zCDP accounting. The yearly mean level is released with separate noise.
- The number of components, the clipping bound, the initialization, and the iteration count are public and not learned from the data.
- Synthetic series are generated from the released parameters only, so generating more of them costs no extra budget.

Limits, stated plainly:

- With about 51 weekly changes per store, privacy at small $\varepsilon$ costs a lot of accuracy. The release uses a single Gaussian ($K = 1$), because mixtures do worse at these budgets. The code supports $K > 1$ for larger datasets.
- Protecting a single underlying weekly sales value is weaker than the stated $\varepsilon$, because that value appears in several rows.
- The clipping bound and $K = 1$ were chosen after looking at results on the private inputs, which is a small leak outside the stated budget.
- Re-running `scripts/dp_release.R` with a different seed spends additional budget on the same data. The ledger in `dp-synth-data/privacy-ledger.csv` records one release per store.

## Repository map

- [`R/dp_gmm.R`](R/dp_gmm.R): DP-EM, zCDP accounting, level release, synthetic generation
- [`tests/testthat/test-dp_gmm.R`](tests/testthat/test-dp_gmm.R): sensitivity, noise calibration, accounting, and recovery tests
- `scripts/dp_evaluate.R`: accuracy versus budget study; writes `dp-results/`
- `scripts/dp_release.R`: the single private release; writes `dp-synth-data/`
- `dp-synth-data/`: private synthetic series, released parameters, and the privacy ledger
- `dp-results/`: evaluation tables. These are computed on the private inputs and are not themselves differentially private.
- `analysis.qmd`, `analysis.html`: original non-private analysis, kept as a baseline
- `raw-csv/`, `raw-data.xlsx`, `recon/`, `synth-data/`: inputs and non-private outputs of the original analysis
- `renamer.sh`: filename normalization helper

## Known problems in the inputs and the original analysis

- The original analysis shifted every store's reconstructed sales to Store 1's yearly mean. The private pipeline uses each store's own mean.
- For Stores 3, 4, and 5, the reconstructed standard deviation matches a different store's row in the workbook summary sheet. The sheets or the labels are probably misassigned; this is flagged in the ledger and needs checking against the source.
- The first-pass inputs labeled `raw` are not original transaction records. They are project inputs, including first-pass synthetic or summarized data.
- The original analysis samples weekly changes independently, which makes simulated years spread out two to three times more than the real data. Real weekly changes are negatively autocorrelated.

## Run it

Install R and the packages `mclust`, `testthat`, `ggplot2`, and `knitr`, plus Quarto to render the report.

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
Rscript scripts/dp_evaluate.R
Rscript scripts/dp_release.R
quarto render dp-analysis.qmd
```

The original analysis additionally needs `tidyverse`, `readxl`, `janitor`, and `stringr`:

```bash
quarto render analysis.qmd
```

Review the diff before committing regenerated output.

## Responsible reuse

Confirm the provenance, confidentiality requirements, and permission to redistribute every input before reusing this repository. A differentially private release only protects the data it was computed from if the inputs stay private. With real data, `raw-csv/`, `raw-data.xlsx`, `recon/`, and the non-private `synth-data/` must not be published.
