describe("lookup table integrity (canonical Kirby 1999 Table 3 design)", {
  key <- get_lookup_table()

  it("exposes the complete 27 x 7 design", {
    expect_s3_class(key, "data.frame")
    expect_equal(nrow(key), 27L)
    expect_true(all(c("questionid", "magnitude", "kindiff", "k_rank",
                      "ss_amount", "ll_amount", "delay") %in% names(key)))
  })

  it("preserves the canonical k-rank row order (MCQ scoring depends on it)", {
    # score_mcq27()/inn()/prop_ss() rely on this exact order; a data-raw rebuild
    # that re-sorts lookup by questionid would silently corrupt scoring.
    expect_equal(
      key$questionid,
      c(13, 1, 9, 20, 6, 17, 26, 24, 12, 22, 16, 15, 3, 10, 2,
        18, 21, 25, 5, 14, 23, 7, 8, 19, 11, 27, 4)
    )
    expect_true(all(diff(key$k_rank) >= 0))
  })

  it("amounts/delays satisfy Kirby Eq.1 against kindiff", {
    k_eq1 <- (key$ll_amount / key$ss_amount - 1) / key$delay
    expect_true(max(abs(k_eq1 - key$kindiff) / key$kindiff) < 1e-5)
    expect_true(all(key$ll_amount > key$ss_amount))
    expect_true(all(key$ss_amount > 0))
    expect_true(all(key$delay > 0))
  })
})
