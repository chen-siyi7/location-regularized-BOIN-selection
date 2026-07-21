# Location-regularized post-trial dose selection

Reproducibility code for the paper *"Stabilizing maximum tolerated dose selection in small
Phase I trials with a location prior."*

The method adds a smooth penalty against boundary doses to the post-trial dose-selection score of
the Bayesian Optimal Interval (BOIN) design. Trial conduct and the toxicity estimator are
unchanged; only the final selection step is modified. This repository reproduces every table and
figure in the paper.

## Requirements

R (>= 4.0). The analysis uses **base R only** — no packages are required to run `reproduce.R`.
`knitr` and `rmarkdown` are needed only to knit the R Markdown report.

## Quick start

From a terminal, in this folder:

```sh
Rscript reproduce.R          # fast defaults (a few minutes)
Rscript reproduce.R full     # manuscript-scale replicates (10,000; slower)
Rscript reproduce.R test     # ~1-minute smoke test
```

`reproduce.R` prints Tables 2–7 to the console and writes Figures 1–5 as PDF files into the
working directory. Alternatively, knit the report from R:

```r
rmarkdown::render("reproduce.Rmd")   # produces reproduce.html
```

## Files

| File | Contents |
|---|---|
| `boin_core.R` | BOIN trial conduct, isotonic estimator, and both selection rules (standard and location-regularized) |
| `designs.R` | Alternative allocation designs: 3+3, CRM, biased-coin, k-in-a-row |
| `comparators.R` | Logistic-regression (Liao-style) and Bayesian dose-response (Wakayama-style) selection |
| `reproduce.R` | Single script: prints Tables 2–7, writes Figures 1–5 |
| `reproduce.Rmd` | R Markdown version producing an HTML report |

## What it reproduces

| Output | Description |
|---|---|
| Table 2 | Fixed case-study scenarios |
| Table 3 | Out-of-sample tuning and validation |
| Table 4 | Parameter sensitivity by true-MTD group |
| Table 5 | Dose-count, target-rate, and sample-size sweeps |
| Table 6 | Transfer across allocation designs (3+3, CRM, biased-coin, k-in-a-row) |
| Table 7 | Model-based selection comparators |
| Figure 1 | Location-penalty shape |
| Figure 2 | Dose-toxicity scenario curves |
| Figure 3 | Gain by true-MTD position |
| Figure 4 | Structural sweeps |
| Figure 5 | Cross-design transfer |

## Reproducibility

All randomness is seeded, so the analysis is fully reproducible: re-running yields identical
numbers. Setting `rho = 0` recovers standard BOIN selection exactly, and the simulator reproduces
the canonical BOIN operating characteristics. Simulation-based quantities are reproduced within
Monte-Carlo error. Set `QUICK <- FALSE` (or run `Rscript reproduce.R full`) to use the
manuscript-scale replicate counts.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`).

## License

Released under the MIT License (`LICENSE`).
