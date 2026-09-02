## Section 6 numerical checks for InsRiskMRVfinite.tex.
## Base R only. Seed: 20260618.
##
## The finite-time ruin event is evaluated pathwise at every claim epoch. This
## is necessary when premium rates are positive, because a line can recover
## between claim epochs as premium is collected.
##
## Full manuscript run (from the section6-reproducibility directory):
##   Rscript code/NumExp_Sec6.R outputs
##
## Optional positive-integer environment variables for a quick test:
##   SEC6_SEED, SEC6_IND_PATHS, SEC6_GAUSSIAN_PATHS, SEC6_CHUNK_SIZE
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
seed_independent <- seed
seed_gaussian <- seed + 1L

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(x) file.path(out_dir, x)

alpha <- 1.5
mean_severity <- alpha / (alpha - 1)
d <- 3L
m <- 2L
T_horizon <- 1
premium_loadings <- c(0.10, 0.20, 0.50)
u_grid <- c(40, 80, 160, 320, 500)
eps_grid <- c(1e-2, 1e-3, 1e-4)

pareto_rnd <- function(n, shape = alpha) {
  if (n == 0L) return(numeric())
  (1 - runif(n))^(-1 / shape)
}

second_largest_3 <- function(M) {
  pmax(pmin(M[, 1L], M[, 2L]),
       pmin(M[, 1L], M[, 3L]),
       pmin(M[, 2L], M[, 3L]))
}

sorted_uniform_times <- function(n_paths, n_arrivals) {
  U <- matrix(runif(n_paths * n_arrivals, 0, T_horizon),
              nrow = n_paths, ncol = n_arrivals)
  if (n_arrivals == 1L) return(U)
  t(apply(U, 1L, sort))
}

add_monte_carlo_columns <- function(scenarios, counts, n_paths, seed_used) {
  scenarios$exceedances <- counts
  scenarios$mc_probability <- counts / n_paths
  scenarios$mc_se <- sqrt(scenarios$mc_probability *
    (1 - scenarios$mc_probability) / n_paths)
  scenarios$ci_low <- pmax(0, scenarios$mc_probability -
                             1.96 * scenarios$mc_se)
  scenarios$ci_high <- pmin(1, scenarios$mc_probability +
                              1.96 * scenarios$mc_se)
  scenarios$mc_to_tail_approximation <- scenarios$mc_probability /
    scenarios$tail_approximation
  scenarios$n_paths <- n_paths
  scenarios$seed <- seed_used
  scenarios
}

## Allocation vectors and tail approximations.

line_rates <- c(1, 2, 4)
alloc_equal <- rep(1 / d, d)

independent_objective <- function(a) {
  sum(combn(seq_len(d), m, FUN = function(S) {
    prod(line_rates[S] * a[S]^(-alpha))
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

make_base_scenarios <- function(allocation_names, allocation_vectors,
                                capital_function) {
  rows <- list()
  for (allocation_name in allocation_names) {
    a <- allocation_vectors[[allocation_name]]
    rows[[length(rows) + 1L]] <- data.frame(
      kind = "u-grid", allocation = allocation_name, value = u_grid,
      epsilon = NA_real_
    )
    rows[[length(rows) + 1L]] <- data.frame(
      kind = "capital", allocation = allocation_name,
      value = vapply(eps_grid, capital_function, numeric(1), a = a),
      epsilon = eps_grid
    )
  }
  base <- do.call(rbind, rows)
  do.call(rbind, lapply(premium_loadings, function(theta) {
    cbind(premium_loading = theta, base)
  }))
}

ind_allocations <- list(optimal = alloc_ind_opt, equal = alloc_equal)
ind_scenarios <- make_base_scenarios(
  names(ind_allocations), ind_allocations, capital_independent
)
ind_scenarios$K <- ifelse(ind_scenarios$allocation == "optimal",
                          K_ind_opt, K_ind_equal)
ind_scenarios$tail_approximation <- ind_scenarios$K *
  ind_scenarios$value^(-m * alpha)
for (j in seq_len(d)) {
  ind_scenarios[[paste0("premium_p", j)]] <- mean_severity *
    (1 + ind_scenarios$premium_loading) * line_rates[j]
}

gau_allocations <- list(equal = alloc_equal, skewed = alloc_skew)
gaussian_scenarios <- make_base_scenarios(
  names(gau_allocations), gau_allocations, capital_gaussian
)
gaussian_scenarios$K <- ifelse(gaussian_scenarios$allocation == "equal",
                               K_gaussian(alloc_equal),
                               K_gaussian(alloc_skew))
gaussian_scenarios$tail_approximation <- gaussian_scenarios$K /
  r_gaussian(gaussian_scenarios$value)
for (j in seq_len(d)) {
  gaussian_scenarios[[paste0("premium_p", j)]] <- mean_severity *
    (1 + gaussian_scenarios$premium_loading) * lambda_g
}

## Exact pathwise simulation at claim epochs.

make_configurations <- function(allocation_vectors, arrival_rates) {
  configs <- expand.grid(
    premium_loading = premium_loadings,
    allocation = names(allocation_vectors),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  configs$a <- lapply(configs$allocation, function(x) allocation_vectors[[x]])
  configs$p <- lapply(configs$premium_loading, function(theta) {
    (1 + theta) * mean_severity * arrival_rates
  })
  configs
}

update_counts <- function(counts, scenarios, configurations, max_statistics) {
  for (cfg in seq_len(nrow(configurations))) {
    idx <- which(
      scenarios$allocation == configurations$allocation[cfg] &
      abs(scenarios$premium_loading -
            configurations$premium_loading[cfg]) < 1e-12
    )
    counts[idx] <- counts[idx] + vapply(
      scenarios$value[idx],
      function(u) sum(max_statistics[, cfg] > u), numeric(1)
    )
  }
  counts
}

simulate_independent_chunk <- function(nr, configurations) {
  N <- rpois(nr, sum(line_rates) * T_horizon)
  max_statistics <- matrix(-Inf, nr, nrow(configurations))
  for (n in sort(unique(N[N > 0L]))) {
    idx <- which(N == n)
    k <- length(idx)
    times <- sorted_uniform_times(k, n)
    labels <- matrix(
      sample.int(d, k * n, replace = TRUE, prob = line_rates),
      nrow = k, ncol = n
    )
    claims <- matrix(pareto_rnd(k * n), nrow = k, ncol = n)
    L <- matrix(0, k, d)
    local_max <- matrix(-Inf, k, nrow(configurations))
    rows <- seq_len(k)
    for (q in seq_len(n)) {
      cell <- cbind(rows, labels[, q])
      L[cell] <- L[cell] + claims[, q]
      for (cfg in seq_len(nrow(configurations))) {
        net_claims <- L - outer(times[, q], configurations$p[[cfg]])
        standardised <- sweep(net_claims, 2L, configurations$a[[cfg]], "/")
        local_max[, cfg] <- pmax(
          local_max[, cfg], second_largest_3(standardised)
        )
      }
    }
    max_statistics[idx, ] <- local_max
  }
  max_statistics
}

rgauss_pareto <- function(n) {
  if (n == 0L) return(matrix(numeric(), 0L, d))
  G <- matrix(rnorm(n * d), nrow = n) %*% chol(Sigma)
  U <- pmin(pnorm(G), 1 - .Machine$double.eps)
  (1 - U)^(-1 / alpha)
}

simulate_gaussian_chunk <- function(nr, configurations) {
  N <- rpois(nr, lambda_g * T_horizon)
  max_statistics <- matrix(-Inf, nr, nrow(configurations))
  for (n in sort(unique(N[N > 0L]))) {
    idx <- which(N == n)
    k <- length(idx)
    times <- sorted_uniform_times(k, n)
    L <- matrix(0, k, d)
    local_max <- matrix(-Inf, k, nrow(configurations))
    for (q in seq_len(n)) {
      L <- L + rgauss_pareto(k)
      for (cfg in seq_len(nrow(configurations))) {
        net_claims <- L - outer(times[, q], configurations$p[[cfg]])
        standardised <- sweep(net_claims, 2L, configurations$a[[cfg]], "/")
        local_max[, cfg] <- pmax(
          local_max[, cfg], second_largest_3(standardised)
        )
      }
    }
    max_statistics[idx, ] <- local_max
  }
  max_statistics
}

run_path_simulation <- function(n_paths, chunk_size, scenarios,
                                configurations, simulator, seed_used) {
  set.seed(seed_used)
  counts <- numeric(nrow(scenarios))
  completed <- 0L
  while (completed < n_paths) {
    nr <- min(chunk_size, n_paths - completed)
    max_statistics <- simulator(nr, configurations)
    counts <- update_counts(counts, scenarios, configurations, max_statistics)
    completed <- completed + nr
  }
  add_monte_carlo_columns(scenarios, counts, n_paths, seed_used)
}

n_independent <- positive_integer_env("SEC6_IND_PATHS", 10000000L)
n_gaussian <- positive_integer_env("SEC6_GAUSSIAN_PATHS", 10000000L)
chunk_size <- positive_integer_env("SEC6_CHUNK_SIZE", 100000L)

ind_configurations <- make_configurations(ind_allocations, line_rates)
ind_scenarios <- run_path_simulation(
  n_independent, chunk_size, ind_scenarios, ind_configurations,
  simulate_independent_chunk, seed_independent
)

gau_configurations <- make_configurations(
  gau_allocations, rep(lambda_g, d)
)
gaussian_scenarios <- run_path_simulation(
  n_gaussian, chunk_size, gaussian_scenarios, gau_configurations,
  simulate_gaussian_chunk, seed_gaussian
)

write.csv(subset(ind_scenarios, kind == "u-grid"),
          out_file("section6_independent_ruin.csv"), row.names = FALSE)
write.csv(subset(ind_scenarios, kind == "capital"),
          out_file("section6_independent_capital.csv"), row.names = FALSE)
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

plot_model <- function(scenarios, allocation_names, panel_label) {
  plot_data <- subset(scenarios, kind == "u-grid")
  y_max <- max(1.25, 1.08 * max(plot_data$mc_to_tail_approximation))
  styles <- data.frame(
    allocation = rep(allocation_names, each = length(premium_loadings)),
    premium_loading = rep(premium_loadings, times = length(allocation_names)),
    pch = rep(c(16, 1), each = length(premium_loadings)),
    lty = rep(seq_along(premium_loadings), times = length(allocation_names))
  )
  first <- TRUE
  for (s in seq_len(nrow(styles))) {
    z <- subset(
      plot_data,
      allocation == styles$allocation[s] &
        abs(premium_loading - styles$premium_loading[s]) < 1e-12
    )
    if (first) {
      plot(z$value, z$mc_to_tail_approximation, type = "b", log = "x",
           pch = styles$pch[s], lty = styles$lty[s], ylim = c(0, y_max),
           xlab = expression("capital " * u),
           ylab = "Monte Carlo / tail approximation")
      first <- FALSE
    } else {
      lines(z$value, z$mc_to_tail_approximation, type = "b",
            pch = styles$pch[s], lty = styles$lty[s])
    }
  }
  abline(h = 1)
  legend(
    "bottomleft",
    paste0(styles$allocation, ", theta=",
           sprintf("%.2f", styles$premium_loading)),
    pch = styles$pch, lty = styles$lty, bty = "n", cex = 0.68
  )
  mtext(panel_label, side = 1, line = 4.0, font = 2, cex = 0.92)
}

pdf(out_file("Sim_fig.pdf"), width = 9.2, height = 4.4, family = "serif")
par(mfrow = c(1, 2), mar = c(5.4, 4.5, 1.0, 0.7),
    mgp = c(2.5, 0.8, 0))
plot_model(ind_scenarios, c("optimal", "equal"),
           "(a) Independent lines")
plot_model(gaussian_scenarios, c("equal", "skewed"),
           "(b) Gaussian-copula claims")
dev.off()

cat("Simulation configuration:\n")
print(c(independent_seed = seed_independent,
        gaussian_seed = seed_gaussian,
        independent_paths = n_independent,
        gaussian_paths = n_gaussian,
        chunk_size = chunk_size))
cat("Premium loadings and rate vectors:\n")
for (theta in premium_loadings) {
  cat("theta =", theta, "; independent p =",
      paste((1 + theta) * mean_severity * line_rates, collapse = ", "),
      "; Gaussian p =",
      paste(rep((1 + theta) * mean_severity * lambda_g, d),
            collapse = ", "), "\n")
}
cat("Independent optimal allocation:\n")
print(alloc_ind_opt)
cat("Independent constants (optimal, equal):\n")
print(c(K_ind_opt, K_ind_equal))
cat("Independent fixed-capital checks:\n")
print(subset(ind_scenarios, kind == "u-grid",
             select = c(premium_loading, allocation, value, exceedances,
                        mc_probability, mc_se,
                        mc_to_tail_approximation)))
cat("Gaussian fixed-capital checks:\n")
print(subset(gaussian_scenarios, kind == "u-grid",
             select = c(premium_loading, allocation, value, exceedances,
                        mc_probability, mc_se,
                        mc_to_tail_approximation)))
