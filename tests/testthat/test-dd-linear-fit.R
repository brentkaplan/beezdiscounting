describe("fit_dd_linear()", {
  toy <- data.frame(
    id = rep(c("a", "b"), each = 3), x = rep(c(1, 7, 30), 2),
    y = c(0.9, 0.5, 0.2, 0.8, 0.4, 0.1)
  )
  it("returns the class and per-subject table equal to the internal helpers", {
    fit <- fit_dd_linear(toy)
    expect_s3_class(fit, "beezdiscounting_linear")
    expect_equal(nrow(fit$subjects), 2L)
    expect_equal(fit$subjects$id, c("a", "b"))
    ya <- log(1 / c(0.9, 0.5, 0.2) - 1) - log(c(1, 7, 30))
    expect_equal(fit$subjects$logk[1], mean(ya))
    expect_equal(fit$subjects$condition, factor(c("(all)", "(all)")))
    expect_null(fit$design$factor)
  })
  it("single condition still fits the RE component (C = 1)", {
    fit <- fit_dd_linear(toy)
    expect_named(fit$re$mu, "(all)"); expect_equal(fit$re$n_delays, 3L)
  })
  it("factors: between-subject levels; mu named by level; unit = subject", {
    d <- rbind(cbind(toy, cond = "ctl"), cbind(toy, cond = "trt"))
    d$id <- paste0(d$id, "_", d$cond)               # distinct subjects per condition
    d$y[d$cond == "trt"] <- d$y[d$cond == "trt"] * 0.8
    fit <- fit_dd_linear(d, factors = "cond")
    expect_equal(nrow(fit$subjects), 4L)
    expect_named(fit$re$mu, c("ctl", "trt"))
    expect_equal(fit$design$levels, c("ctl", "trt"))
    expect_equal(as.integer(fit$re$n_per_condition), c(2L, 2L))
  })
  it("within-subject frames (same id under two levels) are rejected", {
    d <- rbind(cbind(toy, cond = "ctl"), cbind(toy, cond = "trt"))
    expect_error(fit_dd_linear(d, factors = "cond"), "between-subject")
  })
  it("dd_ip: shared validator clamps y > 1 with a warning; clamp default counts 129 boundary points", {
    expect_warning(fit <- fit_dd_linear(dd_ip), "[Cc]lamp")   # dd_ip has y up to 1.13
    expect_equal(fit$transform$boundary, "clamp"); expect_equal(fit$transform$eps, 0.005)
    expect_equal(fit$transform$n_boundary, sum(pmin(dd_ip$y, 1) %in% c(0, 1)))
    expect_equal(fit$transform$n_clamped, fit$transform$n_boundary)   # ruling: only exact 0/1 move
    expect_equal(fit$transform$n_clamped, 129L)
    expect_false(is.null(fit$re)); expect_equal(nrow(fit$subjects), 100L)
  })
  it("dd_ip: eps defaults to 1/(2 ll) when ll is given", {
    fit <- suppressWarnings(fit_dd_linear(dd_ip, ll = 100))
    expect_equal(fit$transform$eps, 0.005)
    fit2 <- suppressWarnings(fit_dd_linear(dd_ip, ll = 1000))
    expect_equal(fit2$transform$eps, 5e-4)
  })
  it("dd_ip: drop mode -> unbalanced -> re NULL with a warning; subjects kept", {
    w <- capture_warnings(fit <- fit_dd_linear(dd_ip, boundary = "drop"))
    expect_true(any(grepl("balanced", w)))
    expect_null(fit$re); expect_equal(fit$transform$n_dropped, fit$transform$n_boundary)
  })
  it("dd_ip: error mode aborts", {
    expect_error(suppressWarnings(fit_dd_linear(dd_ip, boundary = "error")), "equal 0 or 1")
  })
  it("units with < 2 usable points are dropped with a warning naming them", {
    d <- rbind(toy, data.frame(id = "c", x = 1, y = 0.5))
    expect_warning(fit <- fit_dd_linear(d), "fewer than 2")
    expect_equal(fit$subjects$id, c("a", "b"))
  })
  it("rejects more than one factor and non-existent columns", {
    d <- cbind(toy, f1 = "u", f2 = "v")
    expect_error(fit_dd_linear(d, factors = c("f1", "f2")), "exactly one")
    expect_error(fit_dd_linear(toy, factors = "zzz"))
  })
  it("percent scale is divided by 100 before transforming", {
    p <- toy; p$y <- p$y * 100
    expect_warning(fp <- fit_dd_linear(p, response_scale = "percent"))
    expect_equal(fp$subjects$logk, fit_dd_linear(toy)$subjects$logk)
  })
  it("dropping the only subject in a condition drops that condition level everywhere", {
    d <- data.frame(
      id = c("a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2", "c1"),
      x = c(1, 7, 1, 7, 1, 7, 1, 7, 1),
      y = c(0.9, 0.5, 0.85, 0.55, 0.3, 0.1, 0.35, 0.15, 0.5),
      cond = c("A", "A", "A", "A", "B", "B", "B", "B", "C")
    )
    expect_warning(fit <- fit_dd_linear(d, factors = "cond"), "c1")
    expect_equal(fit$design$levels, c("A", "B"))
    expect_named(fit$re$mu, c("A", "B"))
    expect_equal(glance(fit)$n_conditions, 2L)
    expect_warning(pw <- anova(fit, pairwise = TRUE), "g-hat = 0")
    expect_equal(nrow(pw), 1L)
  })
})

describe("fit_dd_linear polish: argument validation and bookkeeping", {
  it("rejects conf_level outside (0, 1), non-scalar, or non-numeric", {
    for (bad in list(0, 1, 1.5, -0.1, c(0.9, 0.95), "a", NA_real_)) {
      expect_error(fit_dd_linear(dd_ip, conf_level = bad), "strictly between")
    }
    expect_s3_class(suppressWarnings(fit_dd_linear(dd_ip, conf_level = 0.9)), "beezdiscounting_linear")
  })
  it("rejects a factor column named like a derived column; other user columns are not carried", {
    sim <- simulate_dd_linear(n_subjects = 4, delays = c(7, 30, 90), mu = c(A = -6, B = -5),
                              sigma2 = 1, g = 5, seed = 11)
    for (nm in c("unit", "y_lin", "d_used", "log_jac")) {
      s <- sim
      s[[nm]] <- s$condition
      expect_error(fit_dd_linear(s, factors = nm), "reserved")
    }
    # a factor literally named `condition` is the normal case; a renamed copy works too
    sim$grp <- sim$condition
    expect_s3_class(fit_dd_linear(sim, factors = "condition"), "beezdiscounting_linear")
    fg <- fit_dd_linear(sim, factors = "grp")
    expect_equal(fg$design$levels, c("A", "B"))
    expect_false("unit_user" %in% names(fg$data))
    # non-factor user columns never reach the fit (validation keeps id/x/y + factor)
    d <- dd_ip
    d$y_lin <- 1
    expect_false(any(suppressWarnings(fit_dd_linear(d))$data$y_lin == 1))
  })
  it("transform counts refer to the retained units (post-drop)", {
    d <- data.frame(
      id = rep(c("a", "b", "c"), each = 3), x = rep(c(1, 7, 30), 3),
      y = c(0.9, 0.5, 0.2, 0, 0, 0.3, 0, 0.5, 0.2)
    )
    fit <- suppressWarnings(fit_dd_linear(d, boundary = "drop"))
    expect_equal(as.character(fit$subjects$id), c("a", "c"))
    expect_equal(nrow(fit$data), 6L)
    expect_equal(fit$transform$n_boundary, 1L)
    expect_equal(fit$transform$n_dropped, 1L)
    expect_equal(fit$transform$n_clamped, 0L)
    expect_output(print(fit), "1 point\\(s\\) at 0/1, 0 clamped, 1 dropped")
  })
  it("balanced boundary = 'drop' keeps the random-effects component", {
    sim <- simulate_dd_linear(n_subjects = 6, delays = c(1, 7, 30, 90), mu = -5,
                              sigma2 = 1, g = 5, seed = 12)
    sim$y[sim$x == 1] <- 0
    w <- capture_warnings(fit <- fit_dd_linear(sim, boundary = "drop"))
    expect_length(w, 0L)
    expect_false(is.null(fit$re))
    expect_equal(fit$re$n_delays, 3L)
    expect_equal(fit$transform$n_dropped, 6L)
  })
  it("returned data carries d_used and log_jac = -ln(d_used - d_used^2)", {
    fit <- suppressWarnings(fit_dd_linear(dd_ip))
    expect_true(all(c("d_used", "log_jac") %in% names(fit$data)))
    ok <- is.finite(fit$data$y_lin)
    expect_equal(fit$data$log_jac[ok], -log(fit$data$d_used[ok] - fit$data$d_used[ok]^2))
  })
})
