golden <- readRDS(testthat::test_path("fixtures", "mcq27-golden.rds"))

test_that("27-item paths are byte-identical to pre-refactor goldens", {
  # class surface is part of the contract
  expect_identical(class(golden$default), c("score_mcq27_output", "data.frame"))
  expect_identical(score_mcq27(mcq27), golden$default)
  expect_identical(score_mcq27(mcq27, trans = "ln"), golden$ln)
  expect_identical(
    score_mcq27(mcq27, trans = "log", round = 3),
    golden$log_round3
  )
  all_sir <- data.frame(subjectid = 1L, questionid = 1:27, response = 0)
  expect_identical(score_mcq27(all_sir), golden$all_sir) # exercises overall+magnitude edges
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
})
# NOTE for implementer: if expect_identical fails ONLY on an attribute-ordering
# artifact (verify with waldo::compare), downgrade that assertion to
# expect_equal(...) plus expect_identical(class(...), ...) and record the
# waldo output in your report. Numeric or structural diffs are real failures.

test_that("score_mcq(items = 27) and the wrapper agree", {
  dat <- generate_data_mcq(n_ids = 6, n_items = 27, seed = 42, prop_na = 0)
  expect_identical(score_mcq27(dat), score_mcq(dat, items = 27))
})

test_that("21-item S3 class is score_mcq_output only (27 keeps legacy class)", {
  dat21 <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  expect_identical(
    class(score_mcq(dat21, items = 21)),
    c("score_mcq_output", "data.frame")
  )
  expect_identical(
    class(score_mcq(mcq27, items = 27)),
    c("score_mcq27_output", "data.frame")
  )
})

test_that("21-item: all-LL and all-SIR edge patterns hit the ladder edges", {
  all_ll <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  res <- score_mcq(all_ll, items = 21)
  expect_equal(res$overall_k, 0.0007, tolerance = 1e-6)
  expect_equal(res$small_k, 0.0007, tolerance = 1e-6)
  expect_equal(res$overall_consistency, 1)
  expect_equal(res$overall_proportion, 1)

  all_sir <- data.frame(subjectid = 1, questionid = 1:21, response = 0)
  res2 <- score_mcq(all_sir, items = 21)
  # overall: geometric midpoint of last item k and edge k (0.1333, scorer "All" sheet)
  expect_equal(res2$overall_k, sqrt(0.1310 * 0.1333), tolerance = 1e-6)
  # per magnitude: edge repeats the magnitude's last kindiff
  expect_equal(res2$small_k, 0.1333, tolerance = 1e-6)
  expect_equal(res2$medium_k, 0.1292, tolerance = 1e-6)
  expect_equal(res2$large_k, 0.1310, tolerance = 1e-6)
  expect_equal(res2$overall_proportion, 0)
})

test_that("21-item: clean switch at Kirby's median (between ladder pos 8 and 9)", {
  # SIR (0) on ladder positions 1-8 (qids 4,15,7,20,9,12,8,16), LL (1) elsewhere.
  # Kirby & Marakovic (1996) report the median subject switches here, k ~ 0.007.
  sir_qids <- c(4, 15, 7, 20, 9, 12, 8, 16)
  dat <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = as.integer(!(1:21 %in% sir_qids))
  )
  res <- score_mcq(dat, items = 21)
  # tolerance widened to 1e-4 on these four: at round = 6 (default) decimal
  # places, a value of this magnitude (~0.007) carries an inherent rounding
  # error of up to ~5e-7 absolute, i.e. ~7e-5 relative -- larger than a 1e-6
  # relative tolerance can absorb. Verified the unrounded value (round = 12)
  # matches the closed form to 9 decimal places, so this is round-off, not a
  # computation error. See task-2-report.md for the full comparison.
  expect_equal(res$overall_k, sqrt(0.0055 * 0.0083), tolerance = 1e-4) # 0.0067565
  expect_equal(res$small_k, sqrt(0.0057 * 0.0083), tolerance = 1e-4)
  expect_equal(res$medium_k, sqrt(0.0055 * 0.0089), tolerance = 1e-4)
  expect_equal(res$large_k, sqrt(0.0031 * 0.0077), tolerance = 1e-4)
  expect_equal(
    res$geomean_k,
    (sqrt(0.0057 * 0.0083) * sqrt(0.0055 * 0.0089) * sqrt(0.0031 * 0.0077))^(1 /
      3),
    tolerance = 1e-6
  )
  expect_equal(res$overall_consistency, 1)
  expect_equal(res$small_consistency, 1)
  expect_equal(res$overall_proportion, 13 / 21, tolerance = 1e-6)
})

test_that("score_mcq validates items argument strictly", {
  expect_error(score_mcq(mcq27, items = 30), "must be 27 or 21")
  expect_error(score_mcq(mcq27, items = 21.9), "must be 27 or 21")
  expect_error(score_mcq(mcq27, items = c(27, 21)), "must be 27 or 21")
})

test_that("score_mcq validates per-subject question coverage and responses", {
  short <- data.frame(subjectid = 1, questionid = 1:20, response = 1)
  expect_error(score_mcq(short, items = 21), "not equal to 21")
  # 21 rows but duplicate id 5 and missing id 6: silent mis-scoring before, error now
  dup <- data.frame(subjectid = 1, questionid = c(1:5, 5, 7:21), response = 1)
  expect_error(score_mcq(dup, items = 21), "not equal to 21")
  unknown <- data.frame(subjectid = 1, questionid = c(1:20, 25), response = 1)
  expect_error(score_mcq(unknown, items = 21), "not equal to 21")
  bad_resp <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = c(2, rep(1, 20))
  )
  expect_error(score_mcq(bad_resp, items = 21), "0, 1, or NA")
  # same guards on the 27-item path
  dup27 <- data.frame(subjectid = 1, questionid = c(1:26, 26), response = 1)
  expect_error(score_mcq27(dup27), "not equal to 27")
})

test_that("score_mcq rejects fractional or non-coercible question ids", {
  # as.integer(as.character()) alone would truncate 1.5 to the valid id 1
  frac21 <- data.frame(subjectid = 1, questionid = c(1.5, 2:21), response = 1)
  expect_error(score_mcq(frac21, items = 21), "not equal to 21")
  bad_str21 <- data.frame(
    subjectid = 1,
    questionid = c("x", as.character(2:21)),
    response = 1
  )
  expect_error(score_mcq(bad_str21, items = 21), "not equal to 21")

  frac27 <- data.frame(subjectid = 1, questionid = c(1.5, 2:27), response = 1)
  expect_error(score_mcq27(frac27), "not equal to 27")
  bad_str27 <- data.frame(
    subjectid = 1,
    questionid = c("x", as.character(2:27)),
    response = 1
  )
  expect_error(score_mcq27(bad_str27), "not equal to 27")
})

test_that("score_mcq normalizes character/factor 0/1 responses", {
  dat_num <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = c(0, rep(1, 20))
  )
  dat_chr <- dat_num
  dat_chr$response <- as.character(dat_chr$response)
  dat_fac <- dat_num
  dat_fac$response <- factor(dat_fac$response, levels = c(0, 1))

  res_num <- score_mcq(dat_num, items = 21)
  res_chr <- score_mcq(dat_chr, items = 21)
  res_fac <- score_mcq(dat_fac, items = 21)
  expect_equal(res_chr$overall_k, res_num$overall_k)
  expect_equal(res_fac$overall_k, res_num$overall_k)

  dat_yes <- dat_num
  dat_yes$response <- as.character(dat_yes$response)
  dat_yes$response[1] <- "yes"
  expect_error(score_mcq(dat_yes, items = 21), "0, 1, or NA")

  dat_two <- dat_num
  dat_two$response[1] <- 2
  expect_error(score_mcq(dat_two, items = 21), "0, 1, or NA")
})

test_that("trans, round, and return_data work for 21 items", {
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  res <- score_mcq(dat, items = 21, trans = "log")
  expect_true("log10_overall_k" %in% names(res))
  expect_equal(res$log10_overall_k, log10(0.0007), tolerance = 1e-5)
  out <- score_mcq(dat, items = 21, return_data = TRUE)
  expect_named(out, c("results", "data"))
  expect_equal(nrow(out$data), 21L)
})

test_that("bundled mcq21 dataset scores to its documented values", {
  expect_equal(nrow(mcq21), 42L)
  res <- score_mcq(mcq21, items = 21)
  expect_equal(
    res$overall_k[res$subjectid == 1],
    sqrt(0.0055 * 0.0083),
    tolerance = 1e-4
  )
  expect_equal(res$overall_k[res$subjectid == 2], 0.0007, tolerance = 1e-6)
  expect_equal(res$overall_consistency, c(1, 1))
})

test_that("score_mcq scores logical responses identically to numeric 0/1 (21-item)", {
  dat_num <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = c(0, rep(1, 20))
  )
  dat_log <- dat_num
  dat_log$response <- as.logical(dat_log$response)

  expect_identical(score_mcq(dat_log, items = 21), score_mcq(dat_num, items = 21))
})

test_that("score_mcq scores logical responses identically to numeric 0/1 (27-item)", {
  dat_log27 <- mcq27
  dat_log27$response <- as.logical(dat_log27$response)

  expect_identical(score_mcq27(dat_log27), score_mcq27(mcq27))
})

test_that("logical NA responses behave like numeric NA responses", {
  dat_num_na <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = c(0, rep(1, 20))
  )
  dat_num_na$response[2] <- NA
  dat_log_na <- dat_num_na
  dat_log_na$response <- as.logical(dat_log_na$response)

  expect_identical(
    score_mcq(dat_log_na, items = 21, impute_method = "none"),
    score_mcq(dat_num_na, items = 21, impute_method = "none")
  )
})
