test_that("score_pdq reproduces Gray et al. (2016) lookups exhaustively (public API)", {
  # The scoring oracle: every one of the 1024 response patterns per block,
  # from the published supplementary lookup tables -- scored END TO END
  # through score_pdq() (validation, merge, block splits, ladder calls),
  # not just through .score_ladder(). Each subject s applies pattern
  # bits(s - 1) within every block. Sequence coding (PDQ_R_Syntax.txt):
  # responses 1 = guaranteed / 2 = risky, item i weighted 2^(i-1),
  # Seq = weighted sum - 1022, so Seq - 1 in binary = the risky (our 1)
  # indicators, least-significant bit = item 1.
  pat <- vapply(
    1:1024,
    function(s) as.integer(intToBits(s - 1))[1:10],
    integer(10)
  ) # 10 x 1024, col = subject
  dat <- data.frame(
    subjectid = rep(1:1024, each = 30),
    questionid = rep(1:30, times = 1024),
    response = as.vector(apply(pat, 2, function(bits) rep(bits, 3)))
  )
  res <- score_pdq(dat, round = 12)
  expect_equal(nrow(res), 1024L)
  expect_identical(res$subjectid, 1:1024)
  for (b in 1:3) {
    lk <- read.table(
      testthat::test_path(
        "fixtures",
        "gray-pdq",
        sprintf("lookup%dPDQ.txt", b)
      ),
      header = TRUE
    )
    expect_equal(nrow(lk), 1024L)
    o <- order(lk[[paste0("Block", b, "Seq")]]) # lookup row for Seq = s
    expect_equal(
      res[[paste0("block", b, "_h")]],
      lk[[paste0("Block", b, "h")]][o],
      tolerance = 1e-6
    )
    # expect_equal() on a vector compares the MEAN relative difference, so a
    # single badly scored pattern can hide inside 1024 correct ones; bound
    # the worst element too (actual max deviation is ~1e-8 per block)
    expect_lt(
      max(abs(
        res[[paste0("block", b, "_h")]] - lk[[paste0("Block", b, "h")]][o]
      )),
      5e-8
    )
    expect_equal(
      res[[paste0("block", b, "_consistency")]],
      lk[[paste0("Block", b, "Cons")]][o],
      tolerance = 1e-9
    )
    expect_equal(
      res[[paste0("block", b, "_proportion")]],
      lk[[paste0("Block", b, "RCR")]][o],
      tolerance = 1e-9
    )
  }
})

test_that("overall_h uses the canonical pooled order at exact-tie items", {
  # q4 and q24 tie exactly at h = 3/4; the canonical pooled order puts q4
  # (position 11) before q24 (position 12). Guaranteed on pooled positions
  # 1-11 (through q4), risky on 12-30 (from q24): a clean pooled switch
  # exactly between the tied pair -> overall_h = sqrt(0.75 * 0.75) = 0.75
  # with consistency 1. A swapped tie order would score this pattern as
  # inconsistent (overall_consistency < 1) -- this pins overall_rank use.
  guar <- c(11, 1, 21, 12, 22, 2, 23, 13, 3, 14, 4)
  dat <- data.frame(
    subjectid = 1,
    questionid = 1:30,
    response = as.integer(!(1:30 %in% guar))
  )
  res <- score_pdq(dat)
  expect_equal(res$overall_h, 0.75, tolerance = 1e-9)
  expect_equal(res$overall_consistency, 1)
})

test_that("score_pdq hand-verified anchors (public API)", {
  all_risky <- data.frame(subjectid = 1, questionid = 1:30, response = 1)
  res <- score_pdq(all_risky)
  expect_identical(class(res), c("score_pdq_output", "data.frame"))
  expect_named(
    res,
    c(
      "subjectid",
      "overall_h",
      "block1_h",
      "block2_h",
      "block3_h",
      "mean_h",
      "geomean_h",
      "overall_consistency",
      "block1_consistency",
      "block2_consistency",
      "block3_consistency",
      "composite_consistency",
      "overall_proportion",
      "block1_proportion",
      "block2_proportion",
      "block3_proportion",
      "impute_method"
    )
  )
  expect_equal(res$block1_h, 0.333333, tolerance = 1e-6)
  expect_equal(res$block2_h, 0.329268, tolerance = 1e-6)
  expect_equal(res$block3_h, 0.333333, tolerance = 1e-6)
  # pooled ladder: q11 (h = 27/82) is the pooled minimum, NOT 1/3
  expect_equal(res$overall_h, 27 / 82, tolerance = 1e-6)
  expect_equal(res$block1_consistency, 1)
  expect_equal(res$overall_consistency, 1)
  expect_equal(res$overall_proportion, 1)

  all_guar <- data.frame(subjectid = 1, questionid = 1:30, response = 0)
  res2 <- score_pdq(all_guar)
  expect_equal(res2$block1_h, 14.647059, tolerance = 1e-6)
  expect_equal(res2$block2_h, 15.166667, tolerance = 1e-6)
  expect_equal(res2$block3_h, 16.166667, tolerance = 1e-6)
  expect_equal(res2$overall_h, 97 / 6, tolerance = 1e-6)
  expect_equal(res2$overall_proportion, 0)

  # clean switch after h rank 5 in every block; also a clean pooled switch
  # (all rank<=5 items pool before all rank>=6 items)
  switch5 <- data.frame(
    subjectid = 1,
    questionid = 1:30,
    response = rep(rep(c(0, 1), each = 5), times = 3)
  )
  res3 <- score_pdq(switch5)
  expect_equal(res3$block1_h, 1.2155706233, tolerance = 1e-6)
  expect_equal(res3$block2_h, 1.2247448714, tolerance = 1e-6)
  expect_equal(res3$block3_h, 1.2339883600, tolerance = 1e-6)
  expect_equal(res3$overall_h, sqrt(1.5), tolerance = 1e-6) # sqrt(67/66 * 99/67)
  expect_equal(res3$mean_h, 1.2247679516, tolerance = 1e-6)
  expect_equal(res3$geomean_h, 1.2247448714, tolerance = 1e-6)
  expect_equal(res3$overall_consistency, 1)
  expect_equal(res3$composite_consistency, 1)
  expect_equal(res3$block2_proportion, 0.5)
})

test_that("score_pdq validates per-subject coverage and responses strictly", {
  short <- data.frame(subjectid = 1, questionid = 1:29, response = 1)
  expect_error(score_pdq(short), "not equal to 30")
  dup <- data.frame(subjectid = 1, questionid = c(1:29, 29), response = 1)
  expect_error(score_pdq(dup), "not equal to 30")
  unknown <- data.frame(subjectid = 1, questionid = c(1:29, 31), response = 1)
  expect_error(score_pdq(unknown), "not equal to 30")
  frac <- data.frame(subjectid = 1, questionid = c(1.5, 2:30), response = 1)
  expect_error(score_pdq(frac), "not equal to 30")
  bad_resp <- data.frame(
    subjectid = 1,
    questionid = 1:30,
    response = c(2, rep(1, 29))
  )
  expect_error(score_pdq(bad_resp), "0, 1, or NA")
})

test_that("score_pdq INN imputes within h-rank groups (one item per block)", {
  dat <- data.frame(subjectid = 1, questionid = 1:30, response = 1)
  dat$response[dat$questionid == 1] <- NA # rank group {1, 11, 21}, both 1
  res <- score_pdq(dat, impute_method = "inn", return_data = TRUE)
  expect_equal(res$data$newresponse[res$data$questionid == 1], 1)
  expect_equal(res$results$block1_h, 0.333333, tolerance = 1e-6)
  # untouched responses round-trip into newresponse
  expect_equal(res$data$newresponse[res$data$questionid != 1], rep(1, 29))
  # uppercase alias dispatches identically
  res_up <- score_pdq(dat, impute_method = "INN", return_data = TRUE)
  expect_equal(res_up$data$newresponse, res$data$newresponse)
  expect_equal(res_up$results$block1_h, res$results$block1_h)

  # disagreeing neighbors leave NA (unless random = TRUE)
  dat2 <- data.frame(subjectid = 1, questionid = 1:30, response = 1)
  dat2$response[dat2$questionid == 11] <- 0
  dat2$response[dat2$questionid == 21] <- NA
  res2 <- score_pdq(dat2, impute_method = "inn", return_data = TRUE)
  expect_true(is.na(res2$data$newresponse[res2$data$questionid == 21]))
  expect_true(is.na(res2$results$block3_h))

  # random = TRUE fills the disagreeing-group NA with a 0/1 draw and
  # labels the method
  set.seed(26)
  res3 <- score_pdq(
    dat2,
    impute_method = "inn",
    random = TRUE,
    return_data = TRUE
  )
  expect_true(res3$data$newresponse[res3$data$questionid == 21] %in% c(0, 1))
  expect_false(is.na(res3$results$block3_h))
  expect_identical(res3$results$impute_method, "INN with random")
})

test_that("score_pdq ggm composites drop NA blocks; none propagates NA", {
  # qid 1 missing -> block 1 (and the pooled overall) score NA
  dat <- data.frame(subjectid = 1, questionid = 1:30, response = 1)
  dat$response[dat$questionid == 1] <- NA
  res_none <- score_pdq(dat, impute_method = "none")
  expect_true(is.na(res_none$block1_h))
  expect_true(is.na(res_none$overall_h))
  expect_true(is.na(res_none$mean_h)) # none: NA propagates
  expect_true(is.na(res_none$geomean_h))
  # ggm: composites average the remaining blocks (27/82 and 1/3)
  res_ggm <- score_pdq(dat, impute_method = "ggm")
  expect_true(is.na(res_ggm$block1_h))
  # anchored at round = 12 so the composites are checked against their exact
  # targets rather than against the default round = 6 quantization (which
  # sits ~1e-6 relative away and would make a 1e-6 tolerance a coin flip)
  res_ggm12 <- score_pdq(dat, impute_method = "ggm", round = 12)
  expect_equal(res_ggm12$mean_h, (27 / 82 + 1 / 3) / 2, tolerance = 1e-9)
  expect_equal(res_ggm12$geomean_h, sqrt(27 / 82 * 1 / 3), tolerance = 1e-9)
  # uppercase alias
  res_up <- score_pdq(dat, impute_method = "GGM")
  expect_equal(res_up$mean_h, res_ggm$mean_h)
})

test_that("trans, round, return_data, and multi-subject input work", {
  dat <- data.frame(
    subjectid = rep(1:2, each = 30),
    questionid = rep(1:30, times = 2),
    response = c(rep(1, 30), rep(0, 30))
  )
  res <- score_pdq(dat)
  expect_equal(nrow(res), 2L)
  logres <- score_pdq(dat, trans = "log")
  expect_true(all(
    c("log10_overall_h", "log10_block1_h", "log10_geomean_h") %in% names(logres)
  ))
  expect_equal(logres$log10_block1_h[1], log10(1 / 3), tolerance = 1e-5)
  lnres <- score_pdq(dat, trans = "ln", round = 3)
  expect_true("ln_mean_h" %in% names(lnres))
  out <- score_pdq(dat, return_data = TRUE)
  expect_named(out, c("results", "data"))
  expect_equal(nrow(out$data), 60L)
  p <- plot(res)
  expect_s3_class(p, "ggplot")
})

test_that("bundled pdq dataset scores to its documented values", {
  expect_equal(nrow(pdq), 60L)
  res <- score_pdq(pdq)
  expect_equal(res$block1_h[res$subjectid == 1], 1.2155706233, tolerance = 1e-6)
  expect_equal(
    res$geomean_h[res$subjectid == 1],
    1.2247448714,
    tolerance = 1e-6
  )
  expect_equal(res$overall_h[res$subjectid == 1], sqrt(1.5), tolerance = 1e-6)
  expect_equal(res$block1_h[res$subjectid == 2], 1 / 3, tolerance = 1e-6)
  expect_equal(res$block2_h[res$subjectid == 2], 27 / 82, tolerance = 1e-6)
  expect_equal(res$overall_h[res$subjectid == 2], 27 / 82, tolerance = 1e-6)
  expect_equal(res$composite_consistency, c(1, 1))
  expect_equal(res$overall_consistency, c(1, 1))
})
