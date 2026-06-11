# Generate the brmsfit fixtures for the dd brms method tests -----------------
# Run from the package root (outside renv if devtools is not in the project
# library):  Rscript data-raw/make-dd-brms-fixtures.R
#
# Each block is guarded with file.exists() so the script only builds MISSING
# fixtures; delete a file (or the whole fixtures/brms dir) to regenerate it.

devtools::load_all(".", quiet = TRUE)

fixture_dir <- file.path("tests", "testthat", "fixtures", "brms")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

path_ip <- file.path(fixture_dir, "fit-mazur-beta.rds")
if (!file.exists(path_ip)) {
  set.seed(31)
  delays <- c(1, 7, 30, 90, 180, 365)
  d <- expand.grid(id = factor(1:6), x = delays)
  k_i <- exp(log(0.02) + rnorm(6, 0, 0.4))
  mu <- 1 / (1 + k_i[d$id] * d$x)
  d$y <- pmin(pmax(mu + rnorm(nrow(d), 0, 0.05), 0), 1)

  message("Fitting mazur/beta fixture...")
  fit <- fit_dd_brms(
    d,
    equation = "mazur", family = "beta",
    chains = 2, iter = 500, warmup = 250,
    cores = 2, seed = 42,
    loo = TRUE, verbose = 1
  )
  saveRDS(fit, path_ip, compress = "xz")
}

path_grp <- file.path(fixture_dir, "fit-mazur-beta-group.rds")
if (!file.exists(path_grp)) {
  message("Fitting mazur/beta + group factor fixture...")
  set.seed(33)
  delays <- c(1, 7, 30, 90, 180, 365)
  d2 <- expand.grid(id = factor(1:8), x = delays)
  d2$group <- factor(ifelse(as.integer(d2$id) <= 4, "ctrl", "treat"))
  k_i2 <- exp(log(0.02) + 0.7 * (d2$group == "treat") + rnorm(8, 0, 0.3)[d2$id])
  mu2 <- 1 / (1 + k_i2 * d2$x)
  d2$y <- pmin(pmax(mu2 + rnorm(nrow(d2), 0, 0.05), 0.01), 0.99)

  fit_grp <- fit_dd_brms(
    d2,
    equation = "mazur", family = "beta",
    factors = "group",
    chains = 2, iter = 500, warmup = 250,
    cores = 2, seed = 42,
    loo = FALSE, verbose = 1
  )
  saveRDS(fit_grp, path_grp, compress = "xz")
}

# Structural choice + group factor (TICKET-048)
path_choice_grp <- file.path(fixture_dir, "fit-choice-mazur-group.rds")
if (!file.exists(path_choice_grp)) {
  message("Fitting choice mazur + group factor fixture...")
  set.seed(35)
  d3 <- expand.grid(
    id = seq_len(8),
    delay = c(1, 7, 30, 90, 180),
    rep = 1:4
  )
  d3$group <- factor(ifelse(d3$id <= 4, "ctrl", "treat"))
  d3$ss_amount <- 50
  d3$ll_amount <- 100
  logk_i3 <- log(0.01) + log(4) * (seq_len(8) > 4) + rnorm(8, 0, 0.3)
  k_i3 <- exp(logk_i3)
  D3 <- 1 / (1 + k_i3[d3$id] * d3$delay)
  eta3 <- 3 * ((d3$ll_amount / d3$ss_amount) * D3 - 1)
  d3$choice <- rbinom(nrow(d3), 1, plogis(eta3))
  d3$id <- factor(d3$id)
  d3$rep <- NULL

  fit_choice_grp <- fit_dd_choice_brms(
    d3,
    equation = "mazur",
    factors = "group",
    chains = 2, iter = 500, warmup = 250,
    cores = 2, seed = 42,
    loo = FALSE, verbose = 1
  )
  saveRDS(fit_choice_grp, path_choice_grp, compress = "xz")
}

saveRDS(
  list(
    brms_version = as.character(packageVersion("brms")),
    created = format(Sys.Date())
  ),
  file.path(fixture_dir, "fixture-meta.rds")
)
for (f in c(path_ip, path_grp, path_choice_grp)) {
  message(sprintf("%s: %.2f MB", basename(f), file.size(f) / 1e6))
}
