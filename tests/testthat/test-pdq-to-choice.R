test_that("pdq_to_choice joins the Madden et al. (2009) design", {
  dat <- data.frame(subjectid = 1, questionid = 1:30, response = 1)
  ch <- pdq_to_choice(dat)
  expect_equal(nrow(ch), 30L)
  expect_named(ch, c("id", "sc_amount", "lu_amount", "prob", "theta", "choice"))
  # spot-checks straight from Table 2
  expect_equal(
    unlist(ch[1, c("sc_amount", "lu_amount", "prob")], use.names = FALSE),
    c(20, 80, 0.10)
  )
  expect_equal(
    unlist(ch[30, c("sc_amount", "lu_amount", "prob")], use.names = FALSE),
    c(40, 60, 0.97)
  )
  expect_equal(ch$theta, (1 - ch$prob) / ch$prob, tolerance = 1e-12)
  expect_true(all(ch$choice == 1))
})

test_that("pdq_to_choice rejects out-of-range ids but keeps ragged/NA input", {
  bad <- data.frame(subjectid = 1, questionid = c(1:29, 31), response = 1)
  expect_error(pdq_to_choice(bad), "Invalid questionid")
  frac <- data.frame(subjectid = 1, questionid = c(1.5, 2:30), response = 1)
  expect_error(pdq_to_choice(frac), "Invalid questionid")
  ragged <- data.frame(
    subjectid = 1,
    questionid = c(2, 5, 19, 23, 30),
    response = c(1, 0, NA, 1, 0)
  )
  ch <- pdq_to_choice(ragged)
  expect_equal(nrow(ch), 5L)
  expect_true(is.na(ch$choice[3]))
  badresp <- data.frame(
    subjectid = 1,
    questionid = 1:30,
    response = c(2, rep(1, 29))
  )
  expect_error(pdq_to_choice(badresp), "must be 0/1")
})

test_that("get_lookup_table keeps back-compat and gains instrument keys", {
  expect_identical(get_lookup_table(), beezdiscounting:::lookup)
  expect_identical(get_lookup_table(21), beezdiscounting:::lookup21) # positional
  expect_identical(get_lookup_table(items = 27), beezdiscounting:::lookup)
  expect_identical(
    get_lookup_table(instrument = "pdq"),
    beezdiscounting:::lookup_pdq
  )
  expect_identical(
    get_lookup_table(instrument = "mcq21"),
    beezdiscounting:::lookup21
  )
  expect_error(get_lookup_table(items = 27, instrument = "pdq"), "not both")
  expect_error(get_lookup_table(instrument = "mcq30"), "must be one of")
})
