## Structural and internal-consistency checks for the Section 6 outputs.
## Run from the repository root with:
##   Rscript code/check_outputs.R outputs

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "outputs"

required <- c(
  "section6_allocation_summary.csv",
  "section6_independent_ruin.csv",
  "section6_independent_capital.csv",
  "section6_gaussian_ruin.csv",
  "section6_gaussian_capital.csv",
  "Sim_fig.pdf"
)
paths <- file.path(out_dir, required)
if (any(!file.exists(paths))) {
  stop("Missing output file(s): ",
       paste(required[!file.exists(paths)], collapse = ", "))
}

read_output <- function(name) {
  read.csv(file.path(out_dir, name), stringsAsFactors = FALSE)
}

allocation <- read_output("section6_allocation_summary.csv")
ind_ruin <- read_output("section6_independent_ruin.csv")
ind_capital <- read_output("section6_independent_capital.csv")
gau_ruin <- read_output("section6_gaussian_ruin.csv")
gau_capital <- read_output("section6_gaussian_capital.csv")

stopifnot(
  nrow(allocation) == 4L,
  nrow(ind_ruin) == 30L,
  nrow(ind_capital) == 18L,
  nrow(gau_ruin) == 30L,
  nrow(gau_capital) == 18L,
  all(abs(rowSums(allocation[c("a1", "a2", "a3")]) - 1) < 1e-12),
  all(abs(allocation$K - c(354.539777918584, 378,
                           11.1644101322913, 13.7832223855448)) < 1e-9)
)

common_columns <- c(
  "premium_loading", "kind", "allocation", "value", "epsilon",
  "premium_p1", "premium_p2", "premium_p3", "exceedances",
  "mc_probability", "mc_se", "ci_low", "ci_high", "K",
  "tail_approximation", "mc_to_tail_approximation", "n_paths", "seed"
)
tables <- list(ind_ruin, ind_capital, gau_ruin, gau_capital)
stopifnot(
  length(unique(unlist(lapply(tables, function(x) x$n_paths)))) == 1L
)
for (x in tables) {
  stopifnot(
    all(common_columns %in% names(x)),
    all(x$mc_probability >= 0 & x$mc_probability <= 1),
    all(x$mc_se >= 0),
    all(x$tail_approximation > 0),
    all(x$premium_loading %in% c(0.1, 0.2, 0.5)),
    all(x$premium_p1 > 0 & x$premium_p2 > 0 & x$premium_p3 > 0),
    all(abs(x$mc_probability - x$exceedances / x$n_paths) < 1e-15),
    all(abs(x$mc_to_tail_approximation -
              x$mc_probability / x$tail_approximation) < 1e-10)
  )
}

for (x in list(ind_ruin, ind_capital)) {
  stopifnot(
    all(abs(x$premium_p1 - 3 * (1 + x$premium_loading)) < 1e-12),
    all(abs(x$premium_p2 - 6 * (1 + x$premium_loading)) < 1e-12),
    all(abs(x$premium_p3 - 12 * (1 + x$premium_loading)) < 1e-12)
  )
}
for (x in list(gau_ruin, gau_capital)) {
  stopifnot(
    all(abs(x$premium_p1 - 3 * (1 + x$premium_loading)) < 1e-12),
    all(abs(x$premium_p2 - x$premium_p1) < 1e-12),
    all(abs(x$premium_p3 - x$premium_p1) < 1e-12)
  )
}

message("All Section 6 output checks passed for: ", normalizePath(out_dir))
