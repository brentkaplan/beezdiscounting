describe(".dd_validate_choice()", {
  good <- data.frame(
    id = c("a", "a", "b", "b"),
    ss_amount = c(40, 55, 31, 14),
    ll_amount = c(65, 75, 85, 25),
    delay = c(25, 61, 14, 19),
    choice = c(0, 1, 1, 0)
  )

  it("returns canonical columns and a coercion_info list on clean data", {
    v <- .dd_validate_choice(good)
    expect_named(v, c("data", "coercion_info"))
    expect_named(v$data, c("id", "ss_amount", "ll_amount", "delay", "choice"))
    expect_equal(nrow(v$data), 4L)
    expect_type(v$coercion_info, "list")
  })

  it("errors on non-binary choice", {
    bad <- good; bad$choice[1] <- 2
    expect_error(.dd_validate_choice(bad), "0/1|binary")
  })

  it("errors on non-positive amounts and non-finite/negative delay", {
    expect_error(.dd_validate_choice(transform(good, ss_amount = c(0, 55, 31, 14))),
                 "ss_amount|positive")
    expect_error(.dd_validate_choice(transform(good, ll_amount = c(65, -1, 85, 25))),
                 "ll_amount|positive")
    expect_error(.dd_validate_choice(transform(good, delay = c(25, Inf, 14, 19))),
                 "delay|finite")
    expect_error(.dd_validate_choice(transform(good, delay = c(25, -1, 14, 19))),
                 "delay|> 0|positive")
  })

  it("warns (not errors) when ll_amount <= ss_amount", {
    odd <- good; odd$ll_amount[1] <- odd$ss_amount[1] - 1
    expect_warning(.dd_validate_choice(odd), "ll_amount|larger")
  })

  it("honors custom column names", {
    nd <- data.frame(subj = "a", ss = 40, ll = 65, d = 25, ch = 1)
    v <- .dd_validate_choice(nd, id_var = "subj", ss_var = "ss", ll_var = "ll",
                             delay_var = "d", choice_var = "ch")
    expect_named(v$data, c("id", "ss_amount", "ll_amount", "delay", "choice"))
  })

  it("errors on non-data-frame input and missing required columns", {
    expect_error(.dd_validate_choice(list()), "must be a data frame")
    expect_error(.dd_validate_choice(good[, -1]), "not found")
  })

  it("retains declared extra_cols verbatim and rejects missing/colliding ones", {
    nd <- data.frame(
      id = c("a", "b"), ss_amount = c(40, 31), ll_amount = c(65, 85),
      delay = c(25, 14), choice = c(0, 1),
      cond = factor(c("C1", "C2")), age = c(30, 41)
    )
    v <- .dd_validate_choice(nd, extra_cols = c("cond", "age"))
    expect_true(all(c("cond", "age") %in% names(v$data)))
    expect_s3_class(v$data$cond, "factor")
    expect_equal(levels(v$data$cond), c("C1", "C2"))
    expect_error(.dd_validate_choice(nd, extra_cols = "nope"), "not found")
  })
})

describe(".dd_choice_prepare_data()", {
  it("0-indexes subject_id aligned to sorted subject_levels and drops NA rows", {
    dat <- data.frame(
      id = c("b", "b", "a", "a"),
      ss_amount = c(40, 55, 31, 14), ll_amount = c(65, 75, 85, 25),
      delay = c(25, 61, 14, 19), choice = c(0, 1, 1, NA)
    )
    prep <- .dd_choice_prepare_data(dat)
    expect_equal(prep$subject_levels, c("a", "b"))
    expect_equal(prep$n_subjects, 2L)
    expect_equal(prep$n_obs, 3L)               # the NA-choice row dropped
    expect_equal(min(prep$subject_id), 0L)
    expect_type(prep$subject_id, "integer")
    expect_equal(length(prep$choice), length(prep$subject_id))
    expect_named(prep$data, c("id", "ss_amount", "ll_amount", "delay", "choice"))
  })

  it("retains extra_cols (droplevels factors) and complete-cases over them too", {
    dat <- data.frame(
      id = c("b", "b", "a", "a"),
      ss_amount = c(40, 55, 31, 14), ll_amount = c(65, 75, 85, 25),
      delay = c(25, 61, 14, 19), choice = c(0, 1, 1, 0),
      cond = factor(c("C1", "C2", "C1", NA), levels = c("C1", "C2", "C3"))
    )
    prep <- .dd_choice_prepare_data(dat, extra_cols = "cond")
    expect_true("cond" %in% names(prep$data))
    expect_equal(prep$n_obs, 3L)                       # NA-cond row dropped
    expect_s3_class(prep$data$cond, "factor")
    expect_equal(levels(prep$data$cond), c("C1", "C2")) # droplevels removed C3
  })
})
