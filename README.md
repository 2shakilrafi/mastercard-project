# Gaussian Mixture Modeling Project

An exploratory University of Arkansas and Mastercard project on reconstructing and generating store-level sales series from summary statistics. The analysis uses R and Gaussian mixture models (GMMs) to create synthetic series whose distributions resemble the supplied aggregates.

## Repository map

- [`analysis.qmd`](analysis.qmd) — primary Quarto analysis
- [`analysis.html`](analysis.html) — rendered report
- `raw-csv/` — worksheet extracts used by the analysis
- `recon/` — reconstructed sales paths
- `synth-data/` — generated synthetic series and path summaries
- `raw-data.xlsx` — source workbook used by the analysis
- `renamer.sh` — filename normalization helper

The files labeled `raw` should not automatically be interpreted as original transaction-level records. The workbook and CSVs used here represent project inputs, including first-pass synthetic or summarized data.

## Reproduce the report

Install R, Quarto, and the packages loaded in `analysis.qmd`, including `tidyverse`, `readxl`, `janitor`, `stringr`, and `mclust`. Then run:

```bash
quarto render analysis.qmd
```

The analysis rewrites derived CSV and report assets. Review the diff before committing regenerated output.

## Responsible reuse

Confirm the provenance, confidentiality requirements, and permission to redistribute every input before reusing this repository. Synthetic data can preserve statistical patterns and may still carry disclosure risk; it should not be treated as anonymous by default.
