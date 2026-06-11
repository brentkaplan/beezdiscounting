# Generate the brmsfit fixture for the dd brms method tests ------------------
# Run from the package root (outside renv if devtools is not in the project
# library):  Rscript data-raw/make-dd-brms-fixtures.R

devtools::load_all(".", quiet = TRUE)

fixture_dir <- file.path("tests", "testthat", "fixtures", "brms")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

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
saveRDS(fit, file.path(fixture_dir, "fit-mazur-beta.rds"), compress = "xz")

message("Fitting mazur/beta + group factor fixture...")
set.seed(33)
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
saveRDS(fit_grp, file.path(fixture_dir, "fit-mazur-beta-group.rds"),
  compress = "xz")

saveRDS(
  list(
    brms_version = as.character(packageVersion("brms")),
    created = format(Sys.Date())
  ),
  file.path(fixture_dir, "fixture-meta.rds")
)
message(sprintf(
  "fixture size: %.2f MB",
  file.size(file.path(fixture_dir, "fit-mazur-beta.rds")) / 1e6
))
