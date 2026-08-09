test_that("prop_sc pools guaranteed-choice proportions by h rank", {
  # subject 1: guaranteed through rank 5 (response 0), risky after;
  # subject 2: all risky
  dat <- data.frame(
    subjectid = rep(1:2, each = 30),
    questionid = rep(1:30, times = 2),
    response = c(rep(rep(c(0, 1), each = 5), times = 3), rep(1, 30))
  )
  res <- prop_sc(dat)
  expect_s3_class(res, "prop_sc_output")
  expect_equal(nrow(res), 10L)
  expect_named(res, c("h_rank", "prop_sc"))
  expect_equal(res$prop_sc[res$h_rank == 1], 0.5) # subj 1 guaranteed, subj 2 risky
  expect_equal(res$prop_sc[res$h_rank == 6], 0) # both risky
  p <- plot(res)
  expect_s3_class(p, "ggplot")
})

test_that("prop_sc warns on missing data and on non-PDQ question ids", {
  dat <- data.frame(
    subjectid = 1,
    questionid = 1:30,
    response = c(NA, rep(1, 29))
  )
  expect_warning(prop_sc(dat), "Missing data")
  bad <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  expect_warning(prop_sc(bad), "do not exactly match")
  # fractional id must not be silently truncated to a valid id
  frac <- data.frame(subjectid = 1, questionid = c(1.5, 2:30), response = 1)
  expect_warning(prop_sc(frac), "do not exactly match")
  # non-coercible id must not be silently dropped from the set comparison
  chr <- data.frame(
    subjectid = 1,
    questionid = c(as.character(1:30), "x"),
    response = 1
  )
  expect_warning(prop_sc(chr), "do not exactly match")
  # extra out-of-range id on top of a complete canonical set
  extra <- data.frame(subjectid = 1, questionid = c(1:30, 31), response = 1)
  expect_warning(prop_sc(extra), "do not exactly match")
})
