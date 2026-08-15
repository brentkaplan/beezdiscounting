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
