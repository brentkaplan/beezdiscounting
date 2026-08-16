# Heavy-file registry for the BEEZ_FULL_TESTS gate (port of beezdemand's
# TICKET-070 / Q14a mechanism).
#
# `.Rbuildignore` excludes `tests`, so R CMD check never runs this suite; the
# only CI job that does is test-coverage.yaml (covr installs from source).
# Under covr's instrumented (-O0 --coverage) TMB build the Monte-Carlo power
# file's thousands of refits push that job past GitHub's 6-hour limit, so the
# heaviest files run only when BEEZ_FULL_TESTS=true:
#   * .github/workflows/full-tests.yaml -- tri-OS, streaming reporter,
#     workflow_dispatch + weekly cron; must be green on every release SHA;
#   * locally: `BEEZ_FULL_TESTS=true Rscript -e 'devtools::test()'`
#     (this repo has no pre-push hook that sets it).
# Ungated runs (test-coverage.yaml, plain devtools::test()) skip them.
#
# Every file listed here calls .skip_unless_full_tests() at its top; the
# registry test in test-full-tests-registry.R keeps the two in sync.
.beez_full_test_files <- c(
  "dd-power"   # Monte Carlo power / find_n: n_sim = 1200/400/400/200/100/60/40/30 TMB refits
  # no-model-fit tests for the same functions live in test-dd-power-fast.R
)

.skip_unless_full_tests <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("BEEZ_FULL_TESTS"), "true"),
    "heavy test file: runs only with BEEZ_FULL_TESTS=true (full-tests.yaml / local full run)"
  )
}
