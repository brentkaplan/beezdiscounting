# ── Delay discounting ──────────────────────────────────────────────────────────

test_that("normalize_dd_response handles numeric recodes", {
  expect_equal(beezdiscounting:::normalize_dd_response("1"), "ss")
  expect_equal(beezdiscounting:::normalize_dd_response("2"), "ll")
  expect_equal(beezdiscounting:::normalize_dd_response(1),   "ss")
  expect_equal(beezdiscounting:::normalize_dd_response(2),   "ll")
})

test_that("normalize_dd_response handles text exports", {
  expect_equal(beezdiscounting:::normalize_dd_response("I'd rather have $100 now"), "ss")
  expect_equal(beezdiscounting:::normalize_dd_response("$150 in 2 months"),         "ll")
})

test_that("score_dd works with text-style example data", {
  result <- score_dd(five.fivetrial_dd)
  item_rows <- result[result$index %in% paste0("I", seq(1, 31, by = 2)), ]
  expect_true(all(item_rows$response %in% c("ss", "ll")))
  expect_true(all(!is.na(item_rows$kval)))
})

test_that("score_dd: numeric recode I17='1' yields response='ss' and SS kval", {
  df <- tibble::tibble(
    ResponseId = "R_num_ss",
    I1 = NA_character_, I3 = NA_character_, I5 = NA_character_,
    I7 = NA_character_, I9 = NA_character_, I11 = NA_character_,
    I13 = NA_character_, I15 = NA_character_,
    I17 = "1",
    I19 = NA_character_, I21 = NA_character_, I23 = NA_character_,
    I25 = NA_character_, I27 = NA_character_, I29 = NA_character_,
    I31 = NA_character_,
    `Attend-SS` = NA_character_, `Attend-LL` = NA_character_
  )
  result <- score_dd(df)
  expect_equal(result$response, "ss")
  expect_equal(result$kval, 0.039551962, tolerance = 1e-6)  # recodess[9]
})

test_that("score_dd: numeric recode I17='2' yields response='ll' and LL kval", {
  df <- tibble::tibble(
    ResponseId = "R_num_ll",
    I1 = NA_character_, I3 = NA_character_, I5 = NA_character_,
    I7 = NA_character_, I9 = NA_character_, I11 = NA_character_,
    I13 = NA_character_, I15 = NA_character_,
    I17 = "2",
    I19 = NA_character_, I21 = NA_character_, I23 = NA_character_,
    I25 = NA_character_, I27 = NA_character_, I29 = NA_character_,
    I31 = NA_character_,
    `Attend-SS` = NA_character_, `Attend-LL` = NA_character_
  )
  result <- score_dd(df)
  expect_equal(result$response, "ll")
  expect_equal(result$kval, 0.023229526, tolerance = 1e-6)  # recodell[9]
})

test_that("score_dd attention check works with numeric recode '1' for Attend-SS", {
  df <- tibble::tibble(
    ResponseId = "R_attn",
    I1 = NA_character_, I3 = NA_character_, I5 = NA_character_,
    I7 = NA_character_, I9 = NA_character_, I11 = NA_character_,
    I13 = NA_character_, I15 = NA_character_, I17 = NA_character_,
    I19 = NA_character_, I21 = NA_character_, I23 = NA_character_,
    I25 = NA_character_, I27 = NA_character_, I29 = NA_character_,
    I31 = NA_character_,
    `Attend-SS` = "1",  # numeric "1" = ss = correct attention response
    `Attend-LL` = NA_character_
  )
  result <- score_dd(df)
  expect_equal(result$attentionflag[result$index == "AttendSS"], "Yes")
})

# ── Probability discounting ────────────────────────────────────────────────────

test_that("normalize_pd_response handles numeric recodes", {
  expect_equal(beezdiscounting:::normalize_pd_response("1"), "sc")
  expect_equal(beezdiscounting:::normalize_pd_response("2"), "lu")
  expect_equal(beezdiscounting:::normalize_pd_response(1),   "sc")
  expect_equal(beezdiscounting:::normalize_pd_response(2),   "lu")
})

test_that("normalize_pd_response handles text exports", {
  expect_equal(beezdiscounting:::normalize_pd_response("$100 for sure"),   "sc")
  expect_equal(beezdiscounting:::normalize_pd_response("50% chance $200"), "lu")
})

test_that("score_pd works with text-style example data", {
  result <- score_pd(five.fivetrial_pd)
  item_rows <- result[result$index %in% paste0("I", seq(1, 31, by = 2)), ]
  expect_true(all(item_rows$response %in% c("sc", "lu")))
  expect_true(all(!is.na(item_rows$hval)))
})

test_that("score_pd: numeric recode I17='1' yields response='sc' and SC hval", {
  df <- tibble::tibble(
    ResponseId = "R_pd_sc",
    I1 = NA_character_, I3 = NA_character_, I5 = NA_character_,
    I7 = NA_character_, I9 = NA_character_, I11 = NA_character_,
    I13 = NA_character_, I15 = NA_character_,
    I17 = "1",
    I19 = NA_character_, I21 = NA_character_, I23 = NA_character_,
    I25 = NA_character_, I27 = NA_character_, I29 = NA_character_,
    I31 = NA_character_,
    `Attend-SS` = NA_character_, `Attend-LL` = NA_character_
  )
  result <- score_pd(df)
  expect_equal(result$response, "sc")
  expect_equal(result$hval, 0.850963, tolerance = 1e-4)  # recodesc[9]
})

test_that("score_pd: numeric recode I17='2' yields response='lu' and LU hval", {
  df <- tibble::tibble(
    ResponseId = "R_pd_lu",
    I1 = NA_character_, I3 = NA_character_, I5 = NA_character_,
    I7 = NA_character_, I9 = NA_character_, I11 = NA_character_,
    I13 = NA_character_, I15 = NA_character_,
    I17 = "2",
    I19 = NA_character_, I21 = NA_character_, I23 = NA_character_,
    I25 = NA_character_, I27 = NA_character_, I29 = NA_character_,
    I31 = NA_character_,
    `Attend-SS` = NA_character_, `Attend-LL` = NA_character_
  )
  result <- score_pd(df)
  expect_equal(result$response, "lu")
  expect_equal(result$hval, 0.624436, tolerance = 1e-4)  # recodelu[9]
})
