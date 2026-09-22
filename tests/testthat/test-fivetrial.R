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

# ── Audit 2026-09-06 F-BZ2-3 / F-BZ2-2: item-aware, closed-vocabulary normalising ──

dd_template_choices <- function() {
  utils::read.csv(test_path("fixtures", "fivetrial", "dd_template_choices.csv"),
                  stringsAsFactors = FALSE)
}

dd_one_row <- function(...) {
  base <- stats::setNames(
    as.list(rep(NA_character_, 33)),
    c(paste0("I", 1:31), "Attend-SS", "Attend-LL")
  )
  vals <- list(...)
  base[names(vals)] <- vals
  tibble::as_tibble(c(list(ResponseId = "R_x"), base))
}

test_that("normalize_dd_response maps every template option (text and code) by item", {
  ch <- dd_template_choices()
  expect_equal(nrow(ch), 66L)
  expect_identical(
    beezdiscounting:::normalize_dd_response(ch$text, ch$item), ch$choice
  )
  expect_identical(
    beezdiscounting:::normalize_dd_response(as.character(ch$code), ch$item),
    ch$choice
  )
  # score_dd() renames the attention columns to AttendSS / AttendLL
  att <- ch[ch$item %in% c("Attend-SS", "Attend-LL"), ]
  expect_identical(
    beezdiscounting:::normalize_dd_response(
      as.character(att$code), gsub("-", "", att$item)
    ),
    att$choice
  )
})

test_that("Attend-LL numeric codes are reversed relative to every other item", {
  # template: Attend-LL choice 1 = "... in 25 years" (LL), 2 = "... now" (SS)
  expect_identical(
    beezdiscounting:::normalize_dd_response(c("1", "2"), c("Attend-LL", "Attend-LL")),
    c("ll", "ss")
  )
  # choosing the 25-year option (numeric 1) fails the attention check
  res <- score_dd(dd_one_row(I17 = "1", `Attend-LL` = "1"))
  expect_identical(unique(res$attentionflag), "Yes")
  # choosing "now" (numeric 2) passes it
  res <- score_dd(dd_one_row(I17 = "1", `Attend-LL` = "2"))
  expect_identical(unique(res$attentionflag), "No")
})

test_that("score_dd gives identical results on the numeric export of the bundled data", {
  ch <- dd_template_choices()
  txt <- five.fivetrial_dd
  num <- txt
  items <- c(paste0("I", 1:31), "Attend-SS", "Attend-LL")
  for (it in intersect(items, names(num))) {
    v <- num[[it]]
    hit <- !is.na(v)
    if (!any(hit)) next
    key <- match(paste(it, v[hit]), paste(ch$item, ch$text))
    expect_false(anyNA(key), info = it)
    v[hit] <- as.character(ch$code[key])
    num[[it]] <- v
  }
  expect_identical(score_dd(num), score_dd(txt))
  expect_identical(ans_dd(num), ans_dd(txt))
})

test_that("normalize_dd_response errors on unrecognised values and keeps NA", {
  expect_error(beezdiscounting:::normalize_dd_response("garbage", "I1"),
               class = "beezdiscounting_response_error")
  expect_error(beezdiscounting:::normalize_dd_response(c("1", "-99"), c("I1", "I3")),
               "-99")
  expect_error(beezdiscounting:::normalize_dd_response("3", "I1"),
               class = "beezdiscounting_response_error")
  expect_identical(beezdiscounting:::normalize_dd_response(c(NA, "1"), c("I1", "I1")),
                   c(NA_character_, "ss"))
  # blank / whitespace = unanswered
  expect_identical(beezdiscounting:::normalize_dd_response(c("", "  "), c("I1", "I3")),
                   c(NA_character_, NA_character_))
  expect_error(score_dd(dd_one_row(I17 = "garbage")),
               class = "beezdiscounting_response_error")
})

test_that("score_dd drops blank items instead of scoring them LL", {
  res <- score_dd(dd_one_row(I1 = "", I17 = "1"))
  expect_identical(res$index, "I17")
  expect_identical(res$response, "ss")
  # a blank attention item must not delete the terminal response
  res <- score_dd(dd_one_row(I1 = "1", `Attend-SS` = ""))
  expect_identical(res$index, "I1")
  expect_equal(nrow(ans_dd(dd_one_row(I1 = "", I17 = "1"))), 1L)
})

test_that("normalize_pd_response errors on unrecognised values and keeps NA", {
  expect_error(beezdiscounting:::normalize_pd_response("garbage", "I1"),
               class = "beezdiscounting_response_error")
  expect_identical(beezdiscounting:::normalize_pd_response(c(NA, "", "2"), c("I1", "I1", "I1")),
                   c(NA_character_, NA_character_, "lu"))
  # the bundled PD export's option texts are all recognised
  pd <- five.fivetrial_pd
  cols <- intersect(c(paste0("I", 1:31), "Attend-SS", "Attend-LL"), names(pd))
  for (it in cols) {
    v <- pd[[it]]
    out <- beezdiscounting:::normalize_pd_response(v, rep(it, length(v)))
    expect_true(all(is.na(v) == is.na(out)), info = it)
  }
  res <- score_pd(dd_one_row(I1 = "", I17 = "1"))
  expect_identical(res$index, "I17")
})

# ── PD template codes (55_Trial_Discounting_Probability_Template_100.qsf) ──
# Fixture built by data-raw/fivetrial_pd_template_choices.R from the Qualtrics
# template: every I-item and Attend-SS code 1 = "for sure" (sc), 2 = "chance"
# (lu); Attend-LL lists the 1%-chance option first (1 = lu, 2 = sc).

pd_template_choices <- function() {
  utils::read.csv(test_path("fixtures", "fivetrial", "pd_template_choices.csv"),
                  stringsAsFactors = FALSE)
}

test_that("normalize_pd_response maps every template option (text and code) by item", {
  ch <- pd_template_choices()
  expect_equal(nrow(ch), 66L)
  expect_identical(
    beezdiscounting:::normalize_pd_response(ch$text, ch$item), ch$choice
  )
  expect_identical(
    beezdiscounting:::normalize_pd_response(as.character(ch$code), ch$item),
    ch$choice
  )
  # score_pd() renames the attention columns to AttendSS / AttendLL
  att <- ch[ch$item %in% c("Attend-SS", "Attend-LL"), ]
  expect_identical(
    beezdiscounting:::normalize_pd_response(
      as.character(att$code), gsub("-", "", att$item)
    ),
    att$choice
  )
})

test_that("PD Attend-LL numeric codes are reversed relative to every other item", {
  expect_identical(
    beezdiscounting:::normalize_pd_response(c("1", "2"), c("Attend-LL", "Attend-LL")),
    c("lu", "sc")
  )
  # choosing the 1%-chance option (numeric 1) fails the attention check
  res <- score_pd(dd_one_row(I17 = "1", `Attend-LL` = "1"))
  expect_identical(unique(res$attentionflag), "Yes")
  # choosing the certain amount (numeric 2) passes it
  res <- score_pd(dd_one_row(I17 = "1", `Attend-LL` = "2"))
  expect_identical(unique(res$attentionflag), "No")
})

test_that("score_pd gives identical results on the numeric export of the bundled data", {
  ch <- pd_template_choices()
  txt <- five.fivetrial_pd
  num <- txt
  items <- c(paste0("I", 1:31), "Attend-SS", "Attend-LL")
  n_mapped <- 0L
  for (it in intersect(items, names(num))) {
    v <- num[[it]]
    hit <- !is.na(v) & nzchar(trimws(v))
    if (!any(hit)) next
    key <- match(paste(it, v[hit]), paste(ch$item, ch$text))
    expect_false(anyNA(key), info = it)
    v[hit] <- as.character(ch$code[key])
    num[[it]] <- v
    n_mapped <- n_mapped + sum(hit)
  }
  expect_gt(n_mapped, 0L)
  expect_identical(score_pd(num), score_pd(txt))
  expect_identical(ans_pd(num), ans_pd(txt))
})
