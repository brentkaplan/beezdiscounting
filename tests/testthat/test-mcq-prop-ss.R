test_that("prop_ss 27-item output is byte-identical to the pre-refactor golden", {
  golden <- readRDS(testthat::test_path("fixtures", "mcq27-golden.rds"))
  expect_identical(suppressWarnings(prop_ss(mcq27)), golden$prop_ss)
  expect_null(attr(suppressWarnings(prop_ss(mcq27)), "mcq_items"))
})

test_that("prop_ss items = 21 groups over 7 ranks with correct proportions", {
  sir_qids <- c(4, 15, 7, 20, 9, 12, 8, 16) # SIR on ranks 1-3 fully
  dat <- data.frame(
    subjectid = 1,
    questionid = 1:21,
    response = as.integer(!(1:21 %in% sir_qids))
  )
  res <- prop_ss(dat, items = 21)
  expect_equal(nrow(res), 7L)
  expect_equal(res$prop_ss[res$k_rank == 1], 1) # all 3 rank-1 items SIR
  expect_equal(res$prop_ss[res$k_rank == 3], 1) # both rank-3 items SIR
  expect_equal(res$prop_ss[res$k_rank == 4], 0) # all 4 rank-4 items LL
  expect_identical(attr(res, "mcq_items"), 21L)
})

test_that("prop_ss plots for both versions and for legacy attribute-less objects", {
  p27 <- plot(suppressWarnings(prop_ss(mcq27))) # no attribute -> 27 labels
  expect_s3_class(p27, "ggplot")
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  p21 <- plot(prop_ss(dat, items = 21))
  expect_s3_class(p21, "ggplot")
})

test_that("prop_ss warns when observed question ids don't match items", {
  dat21 <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  expect_warning(prop_ss(dat21, items = 27), "match")
})

test_that("plot.prop_ss_output labels remaining ranks correctly when a rank is absent", {
  dat <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
  res <- prop_ss(dat, items = 21)
  reg <- .mcq_registry(21)

  # first rank absent
  res_no1 <- res[res$k_rank != 1, ]
  class(res_no1) <- class(res)
  attr(res_no1, "mcq_items") <- 21L
  p1 <- plot(res_no1)
  labs1 <- ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$x$get_labels()
  expect_equal(labs1, reg$rank_labels[2:7])

  # middle rank (4) absent
  res_no4 <- res[res$k_rank != 4, ]
  class(res_no4) <- class(res)
  attr(res_no4, "mcq_items") <- 21L
  p2 <- plot(res_no4)
  labs2 <- ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels()
  expect_equal(labs2, reg$rank_labels[c(1, 2, 3, 5, 6, 7)])
})
