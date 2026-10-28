describe("generate_data_mcq() RNG hygiene", {
  it("restores the caller's RNG state", {
    skip_on_cran()
    set.seed(999)
    before <- get(".Random.seed", envir = globalenv())
    invisible(generate_data_mcq(n_ids = 2, n_items = 5, seed = 42))
    expect_identical(get(".Random.seed", envir = globalenv()), before)

    set.seed(999)
    expected <- runif(1)
    set.seed(999)
    invisible(generate_data_mcq(n_ids = 2, n_items = 5, seed = 42))
    expect_identical(runif(1), expected)
  })

  it("output for a given seed is unchanged by the RNG-restore fix", {
    dat <- generate_data_mcq(n_ids = 2, n_items = 5, seed = 42)
    # Value computed BEFORE the RNG-restore fix (from git HEAD's
    # generate_data_mcq(), which called set.seed(seed) without saving/
    # restoring the caller's state) -- pinned so the fix cannot alter
    # generated output.
    expect_identical(
      dat$response,
      c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 1L)
    )
  })
})

describe("wide_to_long_mcq()", {
  # Each wide row pivots to one block of `items` rows, so questionid must count
  # column positions per row. Sizing it by the number of DISTINCT subjectids
  # aborted ("`questionid` must be size N or 1") whenever a subjectid was
  # repeated or missing more than once (shinybeez prod, 2026-10-28).
  wide <- function(ids, items = 27) {
    out <- as.data.frame(matrix(rep(0:1, length.out = length(ids) * items), nrow = length(ids)))
    colnames(out) <- as.character(seq_len(items))
    cbind(data.frame(subjectid = ids, stringsAsFactors = FALSE), out)
  }

  it("round-trips long_to_wide_mcq() output unchanged for unique subjects", {
    long <- generate_data_mcq(n_ids = 3, n_items = 27, seed = 1)
    back <- as.data.frame(wide_to_long_mcq(long_to_wide_mcq(long)))
    expect_equal(back$subjectid, long$subjectid)
    expect_equal(back$questionid, long$questionid)
    expect_equal(back$response, long$response)
  })

  it("numbers items by column position when a subjectid is repeated", {
    out <- wide_to_long_mcq(wide(c(1, 2, 1)))
    expect_equal(nrow(out), 3 * 27)
    expect_equal(out$questionid, rep(1:27, times = 3))
    expect_equal(out$subjectid, rep(c(1, 2, 1), each = 27))
  })

  it("numbers items by column position when subjectids are missing", {
    out <- wide_to_long_mcq(wide(c(NA, NA)))
    expect_equal(nrow(out), 2 * 27)
    expect_equal(out$questionid, rep(1:27, times = 2))
  })

  it("errors clearly when items does not match the item columns", {
    expect_error(wide_to_long_mcq(wide(1:2, items = 21)), "21 item columns")
    expect_no_error(wide_to_long_mcq(wide(1:2, items = 21), items = 21))
  })
})
