## data-raw/make-mcq27-golden-fixtures.R
## Captures pre-generalization score_mcq27()/prop_ss()/mcq27_to_choice()
## behavior as golden fixtures. Run ONCE at the pre-refactor tree; the
## generalized code must reproduce these objects byte-for-byte.
devtools::load_all(".", quiet = TRUE)
dat_na <- generate_data_mcq(n_ids = 8, n_items = 27, seed = 99, prop_na = 0.03)
all_sir <- data.frame(subjectid = 1L, questionid = 1:27, response = 0)
golden <- list(
  default = score_mcq27(mcq27),
  ln = score_mcq27(mcq27, trans = "ln"),
  log_round3 = score_mcq27(mcq27, trans = "log", round = 3),
  all_sir = score_mcq27(all_sir),
  na_none = suppressWarnings(score_mcq27(dat_na, impute_method = "none")),
  na_ggm = suppressWarnings(score_mcq27(dat_na, impute_method = "ggm")),
  na_inn_data = suppressWarnings(score_mcq27(
    dat_na,
    impute_method = "inn",
    return_data = TRUE
  )),
  prop_ss = suppressWarnings(prop_ss(mcq27)),
  to_choice = mcq27_to_choice(mcq27)
)
dir.create("tests/testthat/fixtures", showWarnings = FALSE)
saveRDS(golden, "tests/testthat/fixtures/mcq27-golden.rds", version = 2)
message("golden fixtures written: ", paste(names(golden), collapse = ", "))
