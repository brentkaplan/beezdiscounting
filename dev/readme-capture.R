# dev/readme-capture.R
#
# Regenerates the *baked* outputs and the TMB figure shown in README.Rmd's
# new modeling-tier sections. Those README code chunks are `eval = FALSE`
# (so `devtools::build_readme()` stays fast and needs no Stan), and the real
# results are pasted in beneath them from this script's output.
#
# Run from the repo root with system R (renv bypassed):
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript dev/readme-capture.R
#
# Writes:
#   dev/readme-outputs.txt           captured console output (copy into README.Rmd)
#   man/figures/readme-tmb-fit.png   TMB population-fit figure (ragg device)
#
# Keep examples small + seeded so the baked outputs are reproducible. Text is
# captured first and the figure is generated last (via ragg::agg_png, which
# avoids the systemfonts segfault the default Rscript png device can hit).

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ggplot2)
})

rule <- function(title) {
  cat("\n\n========== ", title, " ==========\n\n", sep = "")
}

set.seed(1)

## ---- Fit everything (no plotting yet) ----------------------------------
sim <- simulate_dd_ip(
  n_subjects = 15,
  delays = c(7, 30, 90, 180, 365, 730),
  log_k_pop = log(0.01),
  sigma_u = 0.5,
  phi = 10,
  family = "sltb",
  equation = "mazur",
  seed = 1
)
fit <- fit_dd_tmb(
  sim,
  equation = "mazur",
  family = "sltb",
  random_effects = k ~ 1
)

ctrl <- simulate_dd_choice(n_subjects = 12, log_k_pop = log(0.01), seed = 7)
trt <- simulate_dd_choice(n_subjects = 12, log_k_pop = log(0.04), seed = 8)
trt$id <- paste0("t", trt$id)
ctrl$group <- "ctrl"
trt$group <- "treat"
choice_data <- rbind(ctrl, trt)
fit_choice <- fit_dd_choice(
  choice_data,
  mode = "structural",
  factors = "group",
  verbose = 0
)

mcq_choice <- head(mcq27_to_choice(mcq27))

ok_brms <- requireNamespace("brms", quietly = TRUE)

## ---- Capture text outputs ----------------------------------------------
sink("dev/readme-outputs.txt")

rule("TMB: summary(fit)")
print(summary(fit))
rule("TMB: tidy(fit)")
print(tidy(fit))
rule("TMB: glance(fit)")
print(glance(fit))
rule("TMB: head(ranef(fit))")
print(head(ranef(fit)))

rule("CHOICE: get_dd_param_emms(fit_choice)")
print(get_dd_param_emms(fit_choice))
rule("CHOICE: get_dd_comparisons(fit_choice)$k$contrasts_ratio")
print(get_dd_comparisons(fit_choice)$k$contrasts_ratio)

rule("MCQ: head(mcq27_to_choice(mcq27))")
print(mcq_choice)

rule("BRMS: available")
cat("brms available:", ok_brms, "\n")
if (ok_brms) {
  # verbose = 0 makes fit_dd_brms set silent/refresh itself; chains/iter/cores/
  # seed pass through to brms::brm(). tryCatch so a brms hiccup can't drop the
  # TMB/choice/MCQ outputs already written above.
  tryCatch(
    {
      fit_b <- fit_dd_brms(
        dd_ip,
        equation = "mazur",
        family = "beta",
        chains = 2,
        iter = 1000,
        cores = 2,
        seed = 1,
        verbose = 0
      )
      rule("BRMS: summary(fit_b)")
      print(summary(fit_b))
      rule("BRMS: tidy(fit_b)")
      print(tidy(fit_b))
    },
    error = function(e) cat("BRMS fit error:", conditionMessage(e), "\n")
  )
}

cat("\n\nDONE\n")
sink()

## ---- Figure LAST -------------------------------------------------------
# The README figure is the package's own plot() method. ggplot2 + systemfonts
# can segfault computing text metrics under headless Rscript, so render onto an
# explicit cairo device (its own font metrics) and print() -- which is clean.
png(
  "man/figures/readme-tmb-fit.png",
  width = 7,
  height = 4.2,
  units = "in",
  res = 110,
  type = "cairo"
)
print(plot(fit, type = "individual"))
dev.off()
cat("Figure written.\n")
