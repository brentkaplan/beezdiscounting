# plot() for beezdiscounting_linear. Fits are closed form and instant, so each
# block builds its own; fit-dependent tests are skip_on_cran for parity with
# test-dd-plot.R.

.lin_plot_fit_f <- function(conf_level = 0.95) {
  sim <- simulate_dd_linear(
    n_subjects = 8, delays = c(7, 30, 90, 365),
    mu = c(A = -6, B = -4.5), sigma2 = 2, g = 8, seed = 1
  )
  fit_dd_linear(sim, factors = "condition", conf_level = conf_level)
}
.lin_plot_fit_1 <- function() suppressWarnings(fit_dd_linear(dd_ip))
.lin_plot_fit_drop <- function() suppressWarnings(fit_dd_linear(dd_ip, boundary = "drop"))
# Factor fit whose random-effects component is unavailable (one row removed ->
# unbalanced after boundary handling).
.lin_plot_fit_dropf <- function() {
  sim <- simulate_dd_linear(
    n_subjects = 6, delays = c(7, 30, 90, 365),
    mu = c(A = -6, B = -4.5), sigma2 = 2, g = 8, seed = 2
  )
  suppressWarnings(fit_dd_linear(sim[-1, ], factors = "condition"))
}
# Rows in the first layer whose geom inherits `geom` (e.g. "GeomPoint").
.lin_layer_data <- function(p, geom) {
  b <- ggplot2::ggplot_build(p)
  i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))[1]
  b$data[[i]]
}
.lin_n_panels <- function(p) {
  nrow(ggplot2::ggplot_build(p)$layout$layout)
}

test_that("plot.beezdiscounting_linear returns a buildable ggplot for every type", {
  skip_on_cran()
  types <- c("population", "individual", "transformed", "parameters", "resid")
  for (fit in list(.lin_plot_fit_f(), .lin_plot_fit_1())) {
    for (ty in types) {
      p <- suppressMessages(plot(fit, type = ty))
      expect_s3_class(p, "ggplot")
      expect_no_error(ggplot2::ggplot_build(p))
    }
  }
  expect_error(plot(.lin_plot_fit_f(), type = "nope"), "should be one of")
})

test_that("population: one curve per condition, or a single 'Population' curve", {
  skip_on_cran()
  fitf <- .lin_plot_fit_f()
  p <- plot(fitf)
  line <- .lin_layer_data(p, "GeomLine")
  expect_equal(length(unique(line$group)), 2L)
  expect_equal(nrow(line), 2L * 200L)
  # curve values are the hyperbola at exp(mu)
  xs <- line$x[line$group == 1]
  expect_equal(line$y[line$group == 1], 1 / (1 + exp(fitf$re$mu[[1]]) * 10^xs), tolerance = 1e-8)
  fit1 <- .lin_plot_fit_1()
  p1 <- plot(fit1, n_points = 50)
  expect_equal(length(unique(.lin_layer_data(p1, "GeomLine")$group)), 1L)
  expect_equal(nrow(.lin_layer_data(p1, "GeomLine")), 50L)
  expect_s3_class(plot(fit1, x_trans = "linear"), "ggplot")
})

test_that("population: observed points are the raw y and can be hidden", {
  skip_on_cran()
  fit <- .lin_plot_fit_1()
  pts <- .lin_layer_data(plot(fit), "GeomPoint")
  expect_equal(nrow(pts), nrow(fit$data))
  expect_setequal(round(pts$y, 8), round(fit$data$y, 8))
  p <- plot(fit, show_observed = FALSE)
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))))
})

test_that("population errors when the random-effects component was not fitted", {
  skip_on_cran()
  expect_error(plot(.lin_plot_fit_drop()), "not fitted")
})

test_that("individual: per-subject hyperbolae under the population curve; ids subsets", {
  skip_on_cran()
  fitf <- .lin_plot_fit_f()
  p <- plot(fitf, type = "individual")
  b <- ggplot2::ggplot_build(p)
  # layers: observed(1), subject lines(2), population lines(3)
  expect_equal(length(unique(b$data[[2]]$group)), nrow(fitf$subjects))
  expect_equal(length(unique(b$data[[3]]$group)), 2L)
  p2 <- plot(fitf, type = "individual", ids = "A_1", n_points = 20)
  ind <- ggplot2::ggplot_build(p2)$data[[2]]
  expect_equal(nrow(ind), 20L)
  k <- fitf$subjects$k[fitf$subjects$id == "A_1"]
  expect_equal(ind$y, 1 / (1 + k * 10^ind$x), tolerance = 1e-8)
  expect_error(plot(fitf, type = "individual", ids = "nope"), "Unknown")
  expect_error(plot(fitf, type = "individual", ids = character(0)), "empty")
})

test_that("individual draws subject curves only when re is NULL", {
  skip_on_cran()
  fit <- .lin_plot_fit_drop()
  p <- plot(fit, type = "individual", n_points = 10)
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[2]]), 10L * nrow(fit$subjects))
  expect_equal(nrow(b$data[[3]]), 0L)
})

test_that("transformed: facets by subject with slope-1 lines at ln k", {
  skip_on_cran()
  fitf <- .lin_plot_fit_f()
  show <- c("A_1", "B_2")
  p <- plot(fitf, type = "transformed", ids = show)
  expect_equal(.lin_n_panels(p), 2L)
  ab <- .lin_layer_data(p, "GeomAbline")
  expect_equal(nrow(ab), 2L)
  expect_true(all(ab$slope == 1))
  expect_setequal(ab$intercept, fitf$subjects$logk[match(show, fitf$subjects$id)])
  pts <- .lin_layer_data(p, "GeomPoint")
  d <- fitf$data[fitf$data$id %in% show, ]
  expect_equal(nrow(pts), nrow(d))
  expect_setequal(round(pts$y, 8), round(log(1 / d$d_used - 1), 8))
  expect_setequal(round(pts$x, 8), round(log(d$x), 8))
  expect_equal(p$labels$x, "ln(delay)")
})

test_that("transformed: default ids cap at 12 subjects with a message; dropped points absent", {
  skip_on_cran()
  fit1 <- .lin_plot_fit_1()
  expect_message(p <- plot(fit1, type = "transformed"), "12")
  expect_equal(.lin_n_panels(p), 12L)
  expect_no_message(plot(.lin_plot_fit_f(), type = "transformed", ids = "A_1"))
  fitd <- .lin_plot_fit_drop()
  ids <- head(as.character(fitd$subjects$id), 3)
  pd <- plot(fitd, type = "transformed", ids = ids)
  d <- fitd$data[fitd$data$id %in% ids, ]
  expect_equal(nrow(.lin_layer_data(pd, "GeomPoint")), sum(is.finite(d$y_lin)))
})

test_that("parameters with a factor: subject k intervals plus geometric-mean k per condition", {
  skip_on_cran()
  fitf <- .lin_plot_fit_f()
  p <- plot(fitf, type = "parameters")
  b <- ggplot2::ggplot_build(p)
  expect_length(p$layers, 2L)
  expect_equal(nrow(b$data[[1]]), nrow(fitf$subjects))
  pop <- b$data[[2]]
  expect_equal(nrow(pop), 2L)
  # y is log10(k): compare exp(mu) and exp(confint) on the log10 scale
  ci <- confint(fitf, level = fitf$conf_level)
  expect_equal(10^pop$y, unname(exp(fitf$re$mu)), tolerance = 1e-8)
  expect_equal(10^pop$ymin, unname(exp(ci[, 1])), tolerance = 1e-8)
  expect_equal(10^pop$ymax, unname(exp(ci[, 2])), tolerance = 1e-8)
  expect_match(p$labels$y, "log scale")
})

test_that("parameters uses the fit's conf_level for the condition intervals", {
  skip_on_cran()
  fit90 <- .lin_plot_fit_f(conf_level = 0.90)
  pop <- ggplot2::ggplot_build(plot(fit90, type = "parameters"))$data[[2]]
  ci90 <- confint(fit90, level = 0.90)
  ci95 <- confint(fit90, level = 0.95)
  expect_equal(10^pop$ymin, unname(exp(ci90[, 1])), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(10^pop$ymin, unname(exp(ci95[, 1])), tolerance = 1e-8)))
})

test_that("parameters with a factor but no re: subjects only, with a message", {
  skip_on_cran()
  fit <- .lin_plot_fit_dropf()
  expect_null(fit$re)
  expect_message(p <- plot(fit, type = "parameters"), "not fitted")
  expect_length(p$layers, 1L)
  expect_s3_class(p, "ggplot")
})

test_that("parameters without a factor is the subject-k caterpillar", {
  skip_on_cran()
  fit1 <- .lin_plot_fit_1()
  p <- plot(fit1, type = "parameters")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[2]]), nrow(fit1$subjects)) # points
  expect_s3_class(p$theme$axis.text.y, "element_blank") # > 30 subjects
  expect_equal(p$labels$y, "Subject (ordered by k)")
})

test_that("resid plots .std_resid against the raw-scale fitted value", {
  skip_on_cran()
  fit <- .lin_plot_fit_f()
  p <- plot(fit, type = "resid")
  pts <- .lin_layer_data(p, "GeomPoint")
  a <- augment(fit)
  expect_equal(nrow(pts), sum(is.finite(a$.std_resid)))
  expect_setequal(round(pts$x, 8), round(a$.fitted[is.finite(a$.std_resid)], 8))
})
