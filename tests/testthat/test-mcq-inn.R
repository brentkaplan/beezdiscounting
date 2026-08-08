test_that("27-item INN full path matches the pre-refactor golden", {
  golden <- readRDS(testthat::test_path("fixtures", "mcq27-golden.rds"))
  dat_na <- generate_data_mcq(
    n_ids = 8,
    n_items = 27,
    seed = 99,
    prop_na = 0.03
  )
  expect_identical(
    suppressWarnings(score_mcq27(
      dat_na,
      impute_method = "inn",
      return_data = TRUE
    )),
    golden$na_inn_data
  )
})

test_that("27-item INN: NA inside an agreeing rank group takes the group value", {
  dat <- mcq27[mcq27$subjectid == 1, ]
  # rank-1 group is questionids 13, 1, 9; make one NA, others agree
  dat$response[dat$questionid %in% 13] <- NA
  agree <- unique(dat$response[dat$questionid %in% c(1, 9)])
  expect_length(agree, 1L)
  out <- score_mcq(dat, items = 27, impute_method = "inn", return_data = TRUE)
  expect_equal(out$data$newresponse[out$data$questionid == 13], agree)
})

test_that("21-item INN imputes within rank groups, including the size-2 rank 3", {
  # rank 3 for 21 items = qids 8 (S) and 16 (M) only
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  dat$response[dat$questionid == 8] <- NA # neighbor qid 16 is 1 -> impute 1
  res <- score_mcq(dat, items = 21, impute_method = "inn", return_data = TRUE)
  expect_equal(res$data$newresponse[res$data$questionid == 8], 1)
  expect_equal(res$results$overall_k, 0.0007, tolerance = 1e-6)
})

test_that("21-item INN with disagreeing neighbors leaves NA unless random", {
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  # rank 4 group: qids 14, 10, 3, 18 -> make them disagree and one NA
  dat$response[dat$questionid == 14] <- 0
  dat$response[dat$questionid == 10] <- NA
  res <- score_mcq(dat, items = 21, impute_method = "inn", return_data = TRUE)
  expect_true(is.na(res$data$newresponse[res$data$questionid == 10]))
})
