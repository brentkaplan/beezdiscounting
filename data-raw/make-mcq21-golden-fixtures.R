## data-raw/make-mcq21-golden-fixtures.R
## Captures pre-Phase-2a 21-item score_mcq()/prop_ss()/mcq_to_choice()
## behavior as golden fixtures. Run ONCE at the pre-refactor tree
## (develop e70b236); the instrument-keyed registry must reproduce these
## objects byte-for-byte. NEVER regenerate against later code.
devtools::load_all(".", quiet = TRUE)
dat21 <- generate_data_mcq(n_ids = 8, n_items = 21, seed = 99, prop_na = 0)
dat21_na <- generate_data_mcq(
  n_ids = 8,
  n_items = 21,
  seed = 99,
  prop_na = 0.03
)
all_sir <- data.frame(subjectid = 1L, questionid = 1:21, response = 0)
golden <- list(
  default = score_mcq(dat21, items = 21),
  ln = score_mcq(dat21, items = 21, trans = "ln"),
  log_round3 = score_mcq(dat21, items = 21, trans = "log", round = 3),
  all_sir = score_mcq(all_sir, items = 21),
  na_none = suppressWarnings(score_mcq(
    dat21_na,
    items = 21,
    impute_method = "none"
  )),
  na_ggm = suppressWarnings(score_mcq(
    dat21_na,
    items = 21,
    impute_method = "ggm"
  )),
  na_inn_data = suppressWarnings(score_mcq(
    dat21_na,
    items = 21,
    impute_method = "inn",
    return_data = TRUE
  )),
  prop_ss = suppressWarnings(prop_ss(dat21, items = 21)),
  to_choice = mcq_to_choice(dat21, items = 21),
  lookup_table = get_lookup_table(items = 21)
)
saveRDS(golden, "tests/testthat/fixtures/mcq21-golden.rds", version = 2)
message(
  "mcq21 golden fixtures written: ",
  paste(names(golden), collapse = ", ")
)
