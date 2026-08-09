test_that(".instrument_registry serves all three instruments and validates keys", {
  r27 <- beezdiscounting:::.instrument_registry("mcq27")
  r21 <- beezdiscounting:::.instrument_registry("mcq21")
  rpd <- beezdiscounting:::.instrument_registry("pdq")
  expect_identical(r27$table, beezdiscounting:::lookup)
  expect_identical(r21$table, beezdiscounting:::lookup21)
  expect_identical(rpd$table, beezdiscounting:::lookup_pdq)
  expect_identical(r27$items, 27L)
  expect_identical(r27$n_mag, 9L)
  expect_identical(r21$items, 21L)
  expect_identical(r21$n_mag, 7L)
  expect_identical(rpd$items, 30L)
  expect_identical(rpd$n_mag, 10L)
  expect_identical(r27$value_col, "kindiff")
  expect_identical(r27$rank_col, "k_rank")
  expect_identical(rpd$value_col, "hindiff")
  expect_identical(rpd$rank_col, "h_rank")
  expect_identical(rpd$param, "h")
  expect_true(is.na(rpd$edge_k))
  expect_null(r27$overall_order_col)
  expect_null(r21$overall_order_col)
  expect_identical(rpd$overall_order_col, "overall_rank")
  expect_error(
    beezdiscounting:::.instrument_registry("mcq30"),
    "must be one of"
  )
  expect_error(
    beezdiscounting:::.instrument_registry(c("pdq", "mcq27")),
    "must be one of"
  )
  expect_error(
    beezdiscounting:::.instrument_registry(NA_character_),
    "must be one of"
  )
})

test_that(".mcq_registry translates and keeps its strict items contract", {
  expect_identical(
    beezdiscounting:::.mcq_registry(27),
    beezdiscounting:::.instrument_registry("mcq27")
  )
  expect_identical(
    beezdiscounting:::.mcq_registry(21),
    beezdiscounting:::.instrument_registry("mcq21")
  )
  expect_error(beezdiscounting:::.mcq_registry(30), "must be 27 or 21")
  expect_error(beezdiscounting:::.mcq_registry(21.9), "must be 27 or 21")
  expect_error(beezdiscounting:::.mcq_registry(c(27, 21)), "must be 27 or 21")
})

test_that("27-item paths are byte-identical to the pre-refactor goldens", {
  golden <- readRDS(testthat::test_path("fixtures", "mcq27-golden.rds"))
  expect_identical(score_mcq27(mcq27), golden$default)
  expect_identical(score_mcq27(mcq27, trans = "ln"), golden$ln)
  expect_identical(
    score_mcq27(mcq27, trans = "log", round = 3),
    golden$log_round3
  )
  all_sir <- data.frame(subjectid = 1L, questionid = 1:27, response = 0)
  expect_identical(score_mcq27(all_sir), golden$all_sir)
  dat_na <- generate_data_mcq(
    n_ids = 8,
    n_items = 27,
    seed = 99,
    prop_na = 0.03
  )
  expect_identical(
    suppressWarnings(score_mcq27(dat_na, impute_method = "none")),
    golden$na_none
  )
  expect_identical(
    suppressWarnings(score_mcq27(dat_na, impute_method = "ggm")),
    golden$na_ggm
  )
  expect_identical(
    suppressWarnings(score_mcq27(
      dat_na,
      impute_method = "inn",
      return_data = TRUE
    )),
    golden$na_inn_data
  )
  expect_identical(suppressWarnings(prop_ss(mcq27)), golden$prop_ss)
  expect_identical(mcq27_to_choice(mcq27), golden$to_choice)
})

test_that("21-item paths are byte-identical to the pre-refactor goldens", {
  golden <- readRDS(testthat::test_path("fixtures", "mcq21-golden.rds"))
  dat21 <- generate_data_mcq(n_ids = 8, n_items = 21, seed = 99, prop_na = 0)
  dat21_na <- generate_data_mcq(
    n_ids = 8,
    n_items = 21,
    seed = 99,
    prop_na = 0.03
  )
  all_sir <- data.frame(subjectid = 1L, questionid = 1:21, response = 0)
  expect_identical(score_mcq(dat21, items = 21), golden$default)
  expect_identical(score_mcq(dat21, items = 21, trans = "ln"), golden$ln)
  expect_identical(
    score_mcq(dat21, items = 21, trans = "log", round = 3),
    golden$log_round3
  )
  expect_identical(score_mcq(all_sir, items = 21), golden$all_sir)
  expect_identical(
    suppressWarnings(score_mcq(dat21_na, items = 21, impute_method = "none")),
    golden$na_none
  )
  expect_identical(
    suppressWarnings(score_mcq(dat21_na, items = 21, impute_method = "ggm")),
    golden$na_ggm
  )
  expect_identical(
    suppressWarnings(score_mcq(
      dat21_na,
      items = 21,
      impute_method = "inn",
      return_data = TRUE
    )),
    golden$na_inn_data
  )
  expect_identical(suppressWarnings(prop_ss(dat21, items = 21)), golden$prop_ss)
  expect_identical(mcq_to_choice(dat21, items = 21), golden$to_choice)
  expect_identical(get_lookup_table(items = 21), golden$lookup_table)
})

test_that(".score_ladder handles the documented edges", {
  sl <- beezdiscounting:::.score_ladder
  vals <- beezdiscounting:::lookup_pdq$hindiff[1:10] # block 1
  all_risky <- sl(rep(1, 10), vals, vals[10])
  expect_equal(all_risky$value, vals[1], tolerance = 1e-9)
  expect_equal(all_risky$consistency, 1)
  expect_equal(all_risky$proportion, 1)
  all_guar <- sl(rep(0, 10), vals, vals[10])
  expect_equal(all_guar$value, vals[10], tolerance = 1e-9)
  expect_equal(all_guar$consistency, 1)
  expect_equal(all_guar$proportion, 0)
  # NA in, NA out across ALL THREE outputs (impute upstream)
  na_out <- sl(c(NA, rep(1, 9)), vals, vals[10])
  expect_true(is.na(na_out$value))
  expect_true(is.na(na_out$consistency))
  expect_true(is.na(na_out$proportion))
})
