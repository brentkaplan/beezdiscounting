test_that("27-item conversion is byte-identical to the pre-refactor golden", {
  golden <- readRDS(testthat::test_path("fixtures", "mcq27-golden.rds"))
  expect_identical(mcq27_to_choice(mcq27), golden$to_choice)
  expect_identical(mcq_to_choice(mcq27, items = 27), golden$to_choice)
})

test_that("mcq_to_choice(items = 21) joins the Kirby & Marakovic design", {
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  ch <- mcq_to_choice(dat, items = 21)
  expect_equal(nrow(ch), 21L)
  expect_named(ch, c("id", "ss_amount", "ll_amount", "delay", "choice"))
  # spot-checks straight from Table 1
  expect_equal(
    unlist(ch[ch$delay == 43, c("ss_amount", "ll_amount")], use.names = FALSE),
    c(34, 35)
  ) # qid 4
  i1 <- which(dat$questionid == 1)
  expect_equal(
    unlist(ch[i1, c("ss_amount", "ll_amount", "delay")], use.names = FALSE),
    c(30, 85, 14)
  ) # qid 1
  expect_true(all(ch$choice == 1))
})

test_that("mcq_to_choice(items = 21) rejects out-of-range ids but keeps ragged/NA input", {
  bad <- data.frame(subjectid = 1, questionid = c(1:20, 25), response = 1)
  expect_error(mcq_to_choice(bad, items = 21), "Invalid questionid")
  # ragged (subject answered only 5 items) and NA choices stay accepted
  ragged <- data.frame(
    subjectid = 1,
    questionid = c(2, 5, 9, 13, 21),
    response = c(1, 0, NA, 1, 0)
  )
  ch <- mcq_to_choice(ragged, items = 21)
  expect_equal(nrow(ch), 5L)
  expect_true(is.na(ch$choice[3]))
})

test_that("get_lookup_table() default is unchanged; items = 21 returns lookup21", {
  expect_identical(get_lookup_table(), beezdiscounting:::lookup)
  expect_identical(get_lookup_table(items = 21), beezdiscounting:::lookup21)
})
