describe("jarvis2019 dataset", {
  # Load directly from the rda so these tests run with devtools::test() (no install needed).
  jarvis2019 <- local({
    e <- new.env(parent = emptyenv())
    load(testthat::test_path("../../data/jarvis2019.rda"), envir = e)
    e$jarvis2019
  })

  it("has the expected names, dims, and column types", {
    expect_identical(names(jarvis2019), c("id", "x", "y"))
    expect_equal(nrow(jarvis2019), 882L)
    expect_equal(ncol(jarvis2019), 3L)
    expect_true(is.factor(jarvis2019$id))
    expect_true(is.numeric(jarvis2019$x))
    expect_true(is.numeric(jarvis2019$y))
  })

  it("has 126 unique subject-records (row-index ids)", {
    expect_equal(nlevels(jarvis2019$id), 126L)
  })

  it("has exactly the 7 expected delays in days", {
    expect_equal(
      sort(unique(jarvis2019$x)),
      c(7, 30, 180, 365, 730, 1460, 2920)
    )
  })

  it("has y values in [0, 1] with boundary observations present", {
    expect_true(min(jarvis2019$y) == 0)
    expect_true(max(jarvis2019$y) == 1)
    expect_true(sum(jarvis2019$y == 0) > 0)
    expect_true(sum(jarvis2019$y == 1) > 0)
  })
})
