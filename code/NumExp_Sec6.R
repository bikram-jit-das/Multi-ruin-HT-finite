## Section 6 numerical checks for InsRiskMRVfinite.tex.
## Base R only.  Seed: 20260618.
##
## Full manuscript run (from the repository root):
##   Rscript code/NumExp_Sec6.R outputs
##
## Optional positive-integer environment variables for a quick test:
##   SEC6_SEED, SEC6_IND_BATCHES, SEC6_IND_REPS_PER_BATCH,
##   SEC6_GAUSSIAN_PATHS, SEC6_CHUNK_SIZE
##
## Outputs:
##   section6_independent_ruin.csv, section6_independent_capital.csv
##   section6_gaussian_ruin.csv, section6_gaussian_capital.csv
##   section6_allocation_summary.csv, Sim_fig.pdf

positive_integer_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.integer(default))
  value <- suppressWarnings(as.numeric(raw))
  if (!is.finite(value) || value < 1 || value != floor(value)) {
    stop(name, " must be a positive integer.")
  }
  as.integer(value)
}

seed <- positive_integer_env("SEC6_SEED", 20260618L)
set.seed(seed)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(x) file.path(out_dir, x)

alpha <- 1.5
d <- 3L
m <- 2L
T_horizon <- 1
u_grid <- c(40, 80, 160, 320)
eps_grid <- c(1e-2, 1e-3, 1e-4)

pareto_rnd <- function(n, shape = alpha) {
  if (n == 0L) return(numeric())
  (1 - runif(n))^(-1 / shape)
}

pareto_survival <- function(x, shape = alpha) {
  ifelse(x <= 1, 1, x^(-shape))
}

at_least_two_probability <- function(p) {
  p[, 1L] * p[, 2L] + p[, 1L] * p[, 3L] + p[, 2L] * p[, 3L] -
    2 * p[, 1L] * p[, 2L] * p[, 3L]
}

## Independent lines: AK conditional estimator.

line_constants <- c(1, 2, 4)
alloc_equal <- rep(1 / d, d)

independent_objective <- function(a) {
  sum(combn(seq_len(d), m, FUN = function(S) {
    prod(line_constants[S] * a[S]^(-alpha))
  }))
}

softmax3 <- function(theta) {
  z <- c(theta, 0)
  z <- exp(z - max(z))
  z / sum(z)
}

fit_ind <- optim(c(0, 0), function(theta) {
  independent_objective(softmax3(theta))
}, method = "BFGS", control = list(reltol = 1e-13, maxit = 2000L))
alloc_ind_opt <- softmax3(fit_ind$par)

K_ind_equal <- independent_objective(alloc_equal)
K_ind_opt <- independent_objective(alloc_ind_opt)

capital_independent <- function(epsilon, a) {
  (independent_objective(a) / epsilon)^(1 / (m * alpha))
}

## Conditional on N=n, use the Asmussen--Kroese estimator
## n Fbar(max(M_{n-1}, x-S_{n-1})).
ak_tail_batch <- function(nrep, thresholds, rateT, shape = alpha) {
  N <- rpois(nrep, rateT)
  estimates <- matrix(0, nrep, length(thresholds))
  for (n in sort(unique(N))) {
    idx <- which(N == n)
    if (n == 0L) next
    if (n == 1L) {
      partial_sum <- numeric(length(idx))
      partial_max <- rep(1, length(idx))
    } else {
      Y <- matrix(pareto_rnd(length(idx) * (n - 1L), shape),
                  nrow = length(idx))
      partial_sum <- rowSums(Y)
      partial_max <- apply(Y, 1L, max)
    }
    for (k in seq_along(thresholds)) {
      cutoff <- pmax(partial_max, thresholds[k] - partial_sum)
      estimates[idx, k] <- n * pareto_survival(cutoff, shape)
    }
  }
  colMeans(estimates)
}

ind_scenarios <- rbind(
  data.frame(kind = "u-grid", allocation = "optimal", value = u_grid,
             epsilon = NA_real_),
  data.frame(kind = "u-grid", allocation = "equal", value = u_grid,
             epsilon = NA_real_),
  data.frame(kind = "capital", allocation = "optimal",
             value = vapply(eps_grid, capital_independent, numeric(1),
                            a = alloc_ind_opt), epsilon = eps_grid),
  data.frame(kind = "capital", allocation = "equal",
             value = vapply(eps_grid, capital_independent, numeric(1),
                            a = alloc_equal), epsilon = eps_grid)
)
ind_scenarios$epsilon[is.na(ind_scenarios$epsilon)] <- NA_real_

n_batch_ind <- positive_integer_env("SEC6_IND_BATCHES", 40L)
n_per_batch_ind <- positive_integer_env("SEC6_IND_REPS_PER_BATCH", 100000L)
ind_batch_estimates <- matrix(NA_real_, n_batch_ind, nrow(ind_scenarios))

for (r in seq_len(n_batch_ind)) {
  for (allocation_name in c("optimal", "equal")) {
    idx <- which(ind_scenarios$allocation == allocation_name)
    a <- if (allocation_name == "optimal") alloc_ind_opt else alloc_equal
    marginal <- sapply(seq_len(d), function(j) {
      ak_tail_batch(n_per_batch_ind,
                    ind_scenarios$value[idx] * a[j],
                    line_constants[j] * T_horizon)
    })
    ind_batch_estimates[r, idx] <- at_least_two_probability(marginal)
  }
}

ind_scenarios$mc_probability <- colMeans(ind_batch_estimates)
ind_scenarios$mc_se <- apply(ind_batch_estimates, 2L, sd) / sqrt(n_batch_ind)
ind_scenarios$ci_low <- pmax(0, ind_scenarios$mc_probability -
                               1.96 * ind_scenarios$mc_se)
ind_scenarios$ci_high <- ind_scenarios$mc_probability +
  1.96 * ind_scenarios$mc_se
ind_scenarios$K <- ifelse(ind_scenarios$allocation == "optimal",
                          K_ind_opt, K_ind_equal)
ind_scenarios$tail_approximation <- ind_scenarios$K *
  ind_scenarios$value^(-m * alpha)
ind_scenarios$mc_to_tail_approximation <- ind_scenarios$mc_probability /
  ind_scenarios$tail_approximation
ind_scenarios$n_batches <- n_batch_ind
ind_scenarios$n_per_batch <- n_per_batch_ind
ind_scenarios$conditional_replications_per_line <-
  n_batch_ind * n_per_batch_ind

write.csv(subset(ind_scenarios, kind == "u-grid"),
          out_file("section6_independent_ruin.csv"), row.names = FALSE)
write.csv(subset(ind_scenarios, kind == "capital"),
          out_file("section6_independent_capital.csv"), row.names = FALSE)

## Gaussian copula: direct compound-Poisson simulation.

rho <- 0.5
lambda_g <- 1
Sigma <- matrix(rho, d, d)
diag(Sigma) <- 1
alloc_skew <- c(0.20, 0.30, 0.50)

gamma_g <- m / (1 + (m - 1) * rho)
h_g <- 1 / (1 + (m - 1) * rho)
q_g <- alpha * h_g
Upsilon_g <- (1 + (m - 1) * rho)^(m - 0.5) /
  ((2 * pi)^(m / 2) * (1 - rho)^((m - 1) / 2))

r_gaussian <- function(u) {
  (2 * pi)^(-gamma_g / 2) *
    (2 * alpha * log(u))^((m - gamma_g) / 2) *
    u^(alpha * gamma_g)
}

K_gaussian <- function(a) {
  lambda_g * T_horizon * Upsilon_g *
    sum(combn(seq_len(d), m, FUN = function(S) prod(a[S]^(-q_g))))
}

capital_gaussian <- function(epsilon, a) {
  target <- K_gaussian(a) / epsilon
  f <- function(u) r_gaussian(u) - target
  upper <- 2
  while (f(upper) < 0) upper <- 2 * upper
  uniroot(f, c(1 + 1e-8, upper), tol = 1e-11)$root
}

gaussian_scenarios <- rbind(
  data.frame(kind = "u-grid", allocation = "equal", value = u_grid,
             epsilon = NA_real_),
  data.frame(kind = "u-grid", allocation = "skewed", value = u_grid,
             epsilon = NA_real_),
  data.frame(kind = "capital", allocation = "equal",
             value = vapply(eps_grid, capital_gaussian, numeric(1),
                            a = alloc_equal), epsilon = eps_grid),
  data.frame(kind = "capital", allocation = "skewed",
             value = vapply(eps_grid, capital_gaussian, numeric(1),
                            a = alloc_skew), epsilon = eps_grid)
)
gaussian_scenarios$epsilon[is.na(gaussian_scenarios$epsilon)] <- NA_real_

rgauss_pareto <- function(n) {
  if (n == 0L) return(matrix(numeric(), 0L, d))
  G <- matrix(rnorm(n * d), nrow = n) %*% chol(Sigma)
  U <- pmin(pnorm(G), 1 - .Machine$double.eps)
  (1 - U)^(-1 / alpha)
}

second_largest_3 <- function(M) {
  rowSums(M) - apply(M, 1L, min) - apply(M, 1L, max)
}

n_gaussian <- positive_integer_env("SEC6_GAUSSIAN_PATHS", 5000000L)
chunk_size <- positive_integer_env("SEC6_CHUNK_SIZE", 100000L)
gau_counts <- numeric(nrow(gaussian_scenarios))
completed <- 0L

while (completed < n_gaussian) {
  nr <- min(chunk_size, n_gaussian - completed)
  N <- rpois(nr, lambda_g * T_horizon)
  L <- matrix(0, nr, d)
  nz <- which(N > 0L)
  if (length(nz)) {
    id <- rep.int(nz, N[nz])
    Z <- rgauss_pareto(length(id))
    sums <- rowsum(Z, id, reorder = FALSE)
    L[as.integer(rownames(sums)), ] <- sums
  }
  for (allocation_name in c("equal", "skewed")) {
    idx <- which(gaussian_scenarios$allocation == allocation_name)
    a <- if (allocation_name == "equal") alloc_equal else alloc_skew
    standardised <- sweep(L, 2L, a, "/")
    statistic <- second_largest_3(standardised)
    gau_counts[idx] <- gau_counts[idx] + vapply(
      gaussian_scenarios$value[idx],
      function(u) sum(statistic > u), numeric(1))
  }
  completed <- completed + nr
}

gaussian_scenarios$exceedances <- gau_counts
gaussian_scenarios$mc_probability <- gau_counts / n_gaussian
gaussian_scenarios$mc_se <- sqrt(gaussian_scenarios$mc_probability *
  (1 - gaussian_scenarios$mc_probability) / n_gaussian)
gaussian_scenarios$ci_low <- pmax(0, gaussian_scenarios$mc_probability -
                                    1.96 * gaussian_scenarios$mc_se)
gaussian_scenarios$ci_high <- gaussian_scenarios$mc_probability +
  1.96 * gaussian_scenarios$mc_se
gaussian_scenarios$K <- ifelse(gaussian_scenarios$allocation == "equal",
                               K_gaussian(alloc_equal),
                               K_gaussian(alloc_skew))
gaussian_scenarios$tail_approximation <- gaussian_scenarios$K /
  r_gaussian(gaussian_scenarios$value)
gaussian_scenarios$mc_to_tail_approximation <-
  gaussian_scenarios$mc_probability / gaussian_scenarios$tail_approximation
gaussian_scenarios$n_paths <- n_gaussian

write.csv(subset(gaussian_scenarios, kind == "u-grid"),
          out_file("section6_gaussian_ruin.csv"), row.names = FALSE)
write.csv(subset(gaussian_scenarios, kind == "capital"),
          out_file("section6_gaussian_capital.csv"), row.names = FALSE)

allocation_summary <- rbind(
  data.frame(model = "Independent lines", allocation = "optimal",
             a1 = alloc_ind_opt[1], a2 = alloc_ind_opt[2],
             a3 = alloc_ind_opt[3], K = K_ind_opt),
  data.frame(model = "Independent lines", allocation = "equal",
             a1 = alloc_equal[1], a2 = alloc_equal[2],
             a3 = alloc_equal[3], K = K_ind_equal),
  data.frame(model = "Gaussian copula", allocation = "equal",
             a1 = alloc_equal[1], a2 = alloc_equal[2],
             a3 = alloc_equal[3], K = K_gaussian(alloc_equal)),
  data.frame(model = "Gaussian copula", allocation = "skewed",
             a1 = alloc_skew[1], a2 = alloc_skew[2],
             a3 = alloc_skew[3], K = K_gaussian(alloc_skew))
)
write.csv(allocation_summary, out_file("section6_allocation_summary.csv"),
          row.names = FALSE)

## Figure used in the manuscript.
pdf(out_file("Sim_fig.pdf"), width = 9.2, height = 4.1, family = "serif")
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 1.0, 0.7),
    mgp = c(2.5, 0.8, 0))

ind_plot <- subset(ind_scenarios, kind == "u-grid")
z <- subset(ind_plot, allocation == "optimal")
plot(z$value, z$mc_to_tail_approximation, type = "b", log = "x", pch = 16,
     ylim = c(0, 5),
     xlab = expression("capital " * u),
     ylab = "Monte Carlo / tail approximation")
z <- subset(ind_plot, allocation == "equal")
lines(z$value, z$mc_to_tail_approximation, type = "b", pch = 1, lty = 2)
abline(h = 1)
legend("topright", c("optimal allocation", "equal allocation"),
       pch = c(16, 1), lty = c(1, 2), bty = "n", cex = 0.83)

gau_plot <- subset(gaussian_scenarios, kind == "u-grid")
z <- subset(gau_plot, allocation == "equal")
plot(z$value, z$mc_to_tail_approximation, type = "b", log = "x", pch = 16,
     ylim = c(0, 5),
     xlab = expression("capital " * u),
     ylab = "Monte Carlo / tail approximation")
z <- subset(gau_plot, allocation == "skewed")
lines(z$value, z$mc_to_tail_approximation, type = "b", pch = 1, lty = 2)
abline(h = 1)
legend("topright", c("equal allocation", "skewed allocation"),
       pch = c(16, 1), lty = c(1, 2), bty = "n", cex = 0.83)
dev.off()

cat("Simulation configuration:\n")
print(c(seed = seed,
        independent_batches = n_batch_ind,
        independent_replications_per_batch = n_per_batch_ind,
        gaussian_paths = n_gaussian,
        gaussian_chunk_size = chunk_size))
cat("Independent optimal allocation:\n")
print(alloc_ind_opt)
cat("Independent constants (optimal, equal):\n")
print(c(K_ind_opt, K_ind_equal))
cat("\nIndependent capital checks:\n")
print(subset(ind_scenarios, kind == "capital",
	             select = c(allocation, epsilon, value, mc_probability,
	                        mc_se, mc_to_tail_approximation)))
cat("\nGaussian constants (equal, skewed):\n")
print(c(K_gaussian(alloc_equal), K_gaussian(alloc_skew)))
cat("\nGaussian capital checks:\n")
print(subset(gaussian_scenarios, kind == "capital",
	             select = c(allocation, epsilon, value, exceedances,
	                        mc_probability, mc_se, mc_to_tail_approximation)))
