describe("mcq27_to_choice()", {

  it("returns the canonical choice frame for the bundled mcq27 data", {
    out <- mcq27_to_choice(mcq27)
    expect_s3_class(out, "tbl_df")
    expect_named(out, c("id", "ss_amount", "ll_amount", "delay", "choice"))
    expect_type(out$id, "character")
    expect_type(out$choice, "double")
    expect_true(all(out$choice %in% c(0, 1)))
    expect_equal(nrow(out), nrow(mcq27))
    expect_true(all(out$ss_amount > 0))
    expect_true(all(out$ll_amount > out$ss_amount))
    expect_true(all(out$delay > 0))
  })

  it("round-trips amounts/delays to the lookup table per questionid (spec 10.6)", {
    out <- mcq27_to_choice(mcq27)
    key <- get_lookup_table()
    i <- match(mcq27$questionid, key$questionid)
    expect_equal(out$ss_amount, key$ss_amount[i])
    expect_equal(out$ll_amount, key$ll_amount[i])
    expect_equal(out$delay, key$delay[i])
    expect_equal(out$choice, as.numeric(mcq27$response))
    expect_equal(out$id, as.character(mcq27$subjectid))
  })

  it("errors on a missing required column, naming it", {
    expect_error(mcq27_to_choice(mcq27[, c("subjectid", "response")]),
                 "questionid|not found")
  })

  it("errors on unknown questionids, listing the offending value", {
    bad <- mcq27
    bad$questionid[1] <- 99
    expect_error(mcq27_to_choice(bad), "99")
  })

  it("errors on non-coercible questionids", {
    bad <- mcq27
    bad$questionid <- as.character(bad$questionid)
    bad$questionid[1] <- "x"
    expect_error(mcq27_to_choice(bad), "x|questionid")
  })

  it("errors on fractional questionids (no silent truncation to a valid item)", {
    num_bad <- mcq27
    num_bad$questionid[1] <- 1.5
    expect_error(mcq27_to_choice(num_bad), "1.5|Invalid|whole")
    chr_bad <- mcq27
    chr_bad$questionid <- as.character(chr_bad$questionid)
    chr_bad$questionid[1] <- "1.5"
    expect_error(mcq27_to_choice(chr_bad), "1.5|Invalid|whole")
  })

  it("accepts character questionids that name valid items", {
    chr <- mcq27
    chr$questionid <- as.character(chr$questionid)
    expect_equal(mcq27_to_choice(chr)$ss_amount, mcq27_to_choice(mcq27)$ss_amount)
  })

  it("coerces factor responses to numeric 0/1", {
    fac <- mcq27
    fac$response <- factor(fac$response, levels = c(0, 1))
    out <- mcq27_to_choice(fac)
    expect_type(out$choice, "double")
    expect_equal(out$choice, as.numeric(mcq27$response))
  })

  it("errors on responses outside {0, 1}", {
    bad <- mcq27
    bad$response[1] <- 2
    expect_error(mcq27_to_choice(bad), "0/1")
  })

  it("errors on non-coercible responses (not silently dropped to NA)", {
    bad <- mcq27
    bad$response <- as.character(bad$response)
    bad$response[1] <- "x"
    expect_error(mcq27_to_choice(bad), "0/1")
  })

  it("honors custom column names", {
    nd <- data.frame(s = "p1", q = 11, r = 1)
    out <- mcq27_to_choice(nd, id_var = "s", question_var = "q", response_var = "r")
    expect_named(out, c("id", "ss_amount", "ll_amount", "delay", "choice"))
    expect_equal(out$id, "p1")
    expect_equal(out$ss_amount, 11)
    expect_equal(out$ll_amount, 30)
    expect_equal(out$delay, 7)
    expect_equal(out$choice, 1)
  })

  it("tolerates ragged input (no 27-item requirement) and preserves NA responses", {
    one <- mcq27[mcq27$subjectid == 1 & mcq27$questionid %in% 1:3, ]
    expect_no_error(mcq27_to_choice(one))
    expect_equal(nrow(mcq27_to_choice(one)), 3L)
    na_in <- one
    na_in$response[1] <- NA
    expect_true(is.na(mcq27_to_choice(na_in)$choice[1]))
  })

  it("errors on non-data-frame input", {
    expect_error(mcq27_to_choice(list()), "data frame")
  })
})

describe("mcq27_to_choice() -> fit_dd_choice() integration", {
  it("produces output that passes the fit_dd_choice() validator", {
    # .dd_validate_choice() is pure R (no TMB), so this runs unconditionally.
    ch <- mcq27_to_choice(mcq27)
    expect_no_error(.dd_validate_choice(as.data.frame(ch)))
  })
})
