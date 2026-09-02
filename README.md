# Reproducibility files for Section 6

This directory reproduces the numerical calculations, four tables, and Figure 1 in Section 6, “Numerical assessment of ruin and allocation formulas.” The computations cover:

1. three independent compound-Poisson claim processes with Pareto(1.5) severities and arrival rates `(1, 2, 4)`; and
2. a common-arrival compound-Poisson process with Pareto(1.5) margins and an equicorrelated Gaussian copula with correlation `0.5`.

Both models are evaluated with positive premium loadings `theta = 0.10`,
`theta = 0.20`, and `theta = 0.50`. The line-specific premium rate is
`p_j = (1 + theta) lambda_j E[Z]`, where `E[Z] = 3`. Thus the independent
model uses premium vectors `(3.3, 6.6, 13.2)`, `(3.6, 7.2, 14.4)`, and
`(4.5, 9, 18)`, while the Gaussian common-arrival model uses
`(3.3, 3.3, 3.3)`, `(3.6, 3.6, 3.6)`, and `(4.5, 4.5, 4.5)`.

The program uses only base R. No contributed R packages are required.
Fixed-capital ruin probabilities are evaluated at
`u = 40, 80, 160, 320, 500`.

## Directory contents

```text
section6-reproducibility/
├── README.md
├── UPLOAD_TO_GITHUB.md
├── CODEBOOK.md
├── .gitignore
├── code/
│   ├── NumExp_Sec6.R
│   └── check_outputs.R
├── environment/
│   └── R-session-info.txt
└── outputs/
    ├── README.md
    ├── section6_allocation_summary.csv
    ├── section6_independent_ruin.csv
    ├── section6_independent_capital.csv
    ├── section6_gaussian_ruin.csv
    ├── section6_gaussian_capital.csv
    └── Sim_fig.pdf
```

The committed files in `outputs/` are the results reported in the manuscript.
They were generated with seeds `20260618` (independent lines) and `20260619`
(Gaussian copula) using R 4.5.2. The four table CSV files also include 95%
Monte Carlo confidence intervals, although Section 6 reports only the
estimates and standard errors.

## Reproduce the manuscript results

Install R, open a terminal in this directory, and run:

```sh
Rscript code/NumExp_Sec6.R outputs
```

The full run uses 10,000,000 direct paths for each model. It overwrites the
six files in `outputs/` with
results from the stated seeds.

The script prints the allocation vectors, leading constants, capital checks, and simulation settings to the terminal. It also regenerates `outputs/Sim_fig.pdf`.

To check that all expected files, rows, columns, probability identities, and allocation constants are present, run:

```sh
Rscript code/check_outputs.R outputs
```

## Quick installation and code check

The following smaller run checks that the code executes, but it is not intended to reproduce the manuscript’s numerical precision:

```sh
SEC6_IND_PATHS=10000 \
SEC6_GAUSSIAN_PATHS=10000 \
SEC6_CHUNK_SIZE=2000 \
Rscript code/NumExp_Sec6.R /tmp/section6-smoke-test
Rscript code/check_outputs.R /tmp/section6-smoke-test
```

On Windows PowerShell, set the same variables with
`$env:SEC6_IND_PATHS="10000"` and analogous commands before running `Rscript`.

## Methods and standard errors

For both models, the script simulates the finite-time ruin event pathwise.
Between claim epochs all three net-claim coordinates decrease because premium
is collected, so ruin can occur only immediately after a claim. At every
claim epoch the program checks whether at least two coordinates satisfy
`L_j(t) - p_j t > u a_j`. This is the exact multidimensional event; it does
not use a marginal inclusion–exclusion approximation.

For a directly simulated sample proportion `psi_hat` based on `n_paths`, the
reported binomial standard error is

```text
sqrt(psi_hat * (1 - psi_hat) / n_paths).
```

The confidence limits in the CSV files are the estimate plus or minus 1.96 standard errors, truncated below at zero.

## Reproducibility notes

- The default random-number seed and all manuscript parameters are set near the top of the R script.
- Optional environment variables are provided only for smoke testing or sensitivity runs; leaving them unset reproduces the manuscript design.
- Run alternative settings into a different output directory so that the committed manuscript outputs are not overwritten.
- Exact last-digit agreement is best checked with the recorded R version. Small platform-dependent numerical differences can occasionally arise from normal-generation or linear-algebra implementations.

See `CODEBOOK.md` for a description of every output column and `UPLOAD_TO_GITHUB.md` for publication instructions.
