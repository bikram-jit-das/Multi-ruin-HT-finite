# Output codebook

## File-to-manuscript mapping

| Output file | Manuscript item |
|---|---|
| `section6_independent_ruin.csv` | Independent-lines ruin-estimate table |
| `section6_independent_capital.csv` | Independent-lines capital-inversion table |
| `section6_gaussian_ruin.csv` | Gaussian-copula ruin-estimate table |
| `section6_gaussian_capital.csv` | Gaussian-copula capital-inversion table |
| `section6_allocation_summary.csv` | Allocation vectors and leading constants used in the text |
| `Sim_fig.pdf` | Figure 1 |

## Common table columns

| Column | Meaning |
|---|---|
| `premium_loading` | Safety loading `theta` in `p_j = (1 + theta) lambda_j E[Z]` |
| `kind` | `u-grid` for a fixed-capital comparison or `capital` for inversion at a target probability |
| `allocation` | Allocation being evaluated |
| `value` | Total capital `u`; for `kind = capital`, this is the approximating capital obtained by inversion |
| `epsilon` | Target ruin probability for a capital-inversion row; missing for fixed-capital rows |
| `premium_p1`, `premium_p2`, `premium_p3` | Premium-rate vector used for the pathwise simulation |
| `mc_probability` | Monte Carlo estimate of the finite-horizon ruin probability |
| `mc_se` | Monte Carlo standard error |
| `ci_low`, `ci_high` | Normal-approximation 95% Monte Carlo confidence limits |
| `K` | Allocation-dependent leading constant in the Section 6 tail formula |
| `tail_approximation` | First-order tail approximation evaluated at `value` |
| `mc_to_tail_approximation` | `mc_probability / tail_approximation` |
| `exceedances` | Number of directly simulated paths satisfying the ruin event |
| `n_paths` | Total number of directly simulated compound-Poisson paths |
| `seed` | Model-specific random-number seed |

For each path, the program checks the simultaneous multidimensional ruin event
at every claim epoch. Because premium makes every net-claim coordinate decrease
between claim epochs, this check is exact over the full horizon.

## Allocation summary columns

`a1`, `a2`, and `a3` are the three reserve-allocation proportions and sum to one. `K` is the corresponding leading ruin constant. The independent model compares the numerically optimal allocation with equal allocation; the Gaussian-copula model compares equal allocation with the deliberately skewed allocation `(0.2, 0.3, 0.5)`.
