# dd brms <-> TMB parity gate (TICKET-040; design section 6.4) -----------------
#
# The gaussian family is the EXACT-parity arm (identical likelihoods):
# posterior medians must sit within 3 TMB SEs of the MLEs. The beta family
# is compared against TMB sltb QUALITATIVELY only (different likelihoods by
# design: beta + squeeze vs scale-location-truncated beta) - logged, not
# asserted. The choice fitter asserts against fit_dd_choice().
#
# Run from the package root (outside renv if devtools is not in the project
# library):  Rscript inst/validation/dd-brms-tmb-parity.R

devtools::load_all(".", quiet = TRUE)
set.seed(20260611)

N_ID <- 20
DELAYS <- c(1, 7, 30, 90, 180, 365)
MCMC <- list(chains = 2, iter = 1500, warmup = 750, cores = 2)

simulate_ip <- function() {
  logk_i <- log(0.02) + rnorm(N_ID, 0, 0.4)
  d <- expand.grid(id = factor(seq_len(N_ID)), x = DELAYS)
  mu <- 1 / (1 + exp(logk_i[d$id]) * d$x)
  d$y <- pmin(pmax(mu + rnorm(nrow(d), 0, 0.06), 0.001), 0.999)
  d
}

compare_pars <- function(label, brms_coefs, tmb_coefs, tmb_se, par_map) {
  rows <- lapply(names(par_map), function(p_brms) {
    p_tmb <- par_map[[p_brms]]
    b <- unname(brms_coefs[names(brms_coefs) == p_brms])[1]
    t_est <- unname(tmb_coefs[names(tmb_coefs) == p_tmb])[1]
    t_se <- unname(tmb_se[names(tmb_se) == p_tmb])[1]
    data.frame(
      model = label,
      parameter = p_brms,
      brms_median = b,
      tmb_mle = t_est,
      tmb_se = t_se,
      abs_diff = abs(b - t_est),
      tol_3se = 3 * t_se,
      pass = abs(b - t_est) < 3 * t_se
    )
  })
  do.call(rbind, rows)
}

results <- list()
d <- simulate_ip()

message("=== mazur / gaussian (exact-parity arm) ===")
tmb_g <- fit_dd_tmb(d, equation = "mazur", family = "gaussian", verbose = 0)
brms_g <- fit_dd_brms(
  d,
  equation = "mazur",
  family = "gaussian",
  chains = MCMC$chains,
  iter = MCMC$iter,
  warmup = MCMC$warmup,
  cores = MCMC$cores,
  seed = 7,
  loo = FALSE,
  verbose = 0
)
results$gaussian <- compare_pars(
  "mazur_gaussian",
  brms_g$model$coefficients,
  tmb_g$model$coefficients,
  tmb_g$model$se,
  list(beta_k = "beta_k")
)
print(results$gaussian, row.names = FALSE, digits = 4)

message("\n=== mazur / beta vs TMB sltb (qualitative; logged only) ===")
tmb_s <- fit_dd_tmb(d, equation = "mazur", family = "sltb", verbose = 0)
brms_b <- fit_dd_brms(
  d,
  equation = "mazur",
  family = "beta",
  chains = MCMC$chains,
  iter = MCMC$iter,
  warmup = MCMC$warmup,
  cores = MCMC$cores,
  seed = 7,
  loo = FALSE,
  verbose = 0
)
message(sprintf(
  "beta logk median = %.4f vs sltb MLE = %.4f (diff %.4f; truth %.4f) - analog, not asserted",
  brms_b$model$coefficients[["beta_k"]],
  tmb_s$model$coefficients[["beta_k"]],
  abs(
    brms_b$model$coefficients[["beta_k"]] - tmb_s$model$coefficients[["beta_k"]]
  ),
  log(0.02)
))

message("\n=== structural choice ===")
sim_choice <- function() {
  logk_i <- log(0.02) + rnorm(N_ID, 0, 0.4)
  d <- expand.grid(
    id = factor(seq_len(N_ID)),
    delay = c(1, 7, 30, 90, 180),
    rep = 1:4
  )
  d$ss_amount <- 50
  d$ll_amount <- 100
  D <- 1 / (1 + exp(logk_i[d$id]) * d$delay)
  d$choice <- rbinom(nrow(d), 1, plogis(3 * (2 * D - 1)))
  d$rep <- NULL
  d
}
dc <- sim_choice()
tmb_c <- fit_dd_choice(dc, mode = "structural", equation = "mazur", verbose = 0)
brms_c <- fit_dd_choice_brms(
  dc,
  equation = "mazur",
  chains = MCMC$chains,
  iter = MCMC$iter,
  warmup = MCMC$warmup,
  cores = MCMC$cores,
  seed = 7,
  loo = FALSE,
  verbose = 0
)
results$choice <- compare_pars(
  "choice_structural",
  brms_c$model$coefficients,
  tmb_c$model$coefficients,
  tmb_c$model$se,
  list(beta_k = "beta_k", log_gamma = "log_gamma")
)
print(results$choice, row.names = FALSE, digits = 4)

all_res <- do.call(rbind, results)
n_fail <- sum(!all_res$pass)
out_path <- file.path("inst", "validation", "dd-brms-tmb-parity-results.csv")
utils::write.csv(all_res, out_path, row.names = FALSE)
message(sprintf(
  "\nParity gate: %d/%d asserted comparisons within 3 TMB SEs (results: %s)",
  sum(all_res$pass),
  nrow(all_res),
  out_path
))
if (n_fail > 0) {
  print(all_res[!all_res$pass, ], row.names = FALSE)
  stop(n_fail, " parity violation(s).", call. = FALSE)
}
message("PARITY GATE PASSED")
