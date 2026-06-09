# Pure-R descriptive-choice helpers: design builder, Cholesky<->Sigma, eta.

describe(".dd_choice_descriptive_design", {
  dat <- data.frame(
    id = rep(c("a", "b"), each = 3),
    ss_amount = c(40, 55, 31, 40, 55, 31),
    ll_amount = c(65, 75, 85, 65, 75, 85),
    delay = c(25, 61, 14, 25, 61, 14),
    choice = c(0, 1, 1, 0, 0, 1)
  )

  it("default design = Young's 2 predictors (no intercept), Z == Zre, q = 2", {
    d <- .dd_choice_descriptive_design(dat, predictors = NULL, random_slopes = TRUE)
    expect_equal(ncol(d$Z), 2L)
    expect_equal(ncol(d$Zre), 2L)
    expect_equal(d$q, 2L)
    expect_equal(unname(d$Z[1, 1]), log(65 / 40), tolerance = 1e-10)
    expect_equal(unname(d$Z[1, 2]), log(25 + 1), tolerance = 1e-10)
    expect_false(any(grepl("Intercept", colnames(d$Z))))
  })

  it("random_slopes = FALSE => q = 0, Zre has 0 columns", {
    d <- .dd_choice_descriptive_design(dat, predictors = NULL, random_slopes = FALSE)
    expect_equal(d$q, 0L)
    expect_equal(ncol(d$Zre), 0L)
    expect_equal(ncol(d$Z), 2L)
  })

  it("predictors override changes Z but NOT Zre (Zre is always Young's 2)", {
    dat2 <- transform(dat, grp = factor(rep(c("x", "y"), 3)))
    d <- .dd_choice_descriptive_design(
      dat2,
      predictors = ~ 0 + log(ll_amount / ss_amount) + log(delay + 1) + grp,
      random_slopes = TRUE
    )
    # ~ 0 + two_continuous + two-level factor => 2 + 2 = 4 cols (no intercept =>
    # both dummy levels appear).  The key invariant is that Zre is unchanged.
    expect_equal(ncol(d$Z), 4L)
    expect_equal(ncol(d$Zre), 2L)
    expect_equal(d$q, 2L)
  })

  it("rejects a rank-deficient fixed design", {
    expect_error(
      .dd_choice_descriptive_design(
        dat,
        predictors = ~ 0 + log(delay + 1) + I(log(delay + 1)),
        random_slopes = TRUE
      ),
      regexp = "rank-deficient|aliased|collinear"
    )
  })
})

describe(".dd_chol_sigma (2-RE)", {
  it("round-trips (sd1, sd2, rho) -> Sigma -> (sd1, sd2, rho)", {
    sd1 <- 0.7; sd2 <- 0.4; rho <- -0.5
    out <- .dd_chol_sigma(log_sd_re = log(c(sd1, sd2)), cor_re = atanh(rho))
    Sigma <- out$Sigma
    expect_equal(sqrt(Sigma[1, 1]), sd1, tolerance = 1e-10)
    expect_equal(sqrt(Sigma[2, 2]), sd2, tolerance = 1e-10)
    expect_equal(Sigma[1, 2] / (sd1 * sd2), rho, tolerance = 1e-10)
    expect_equal(out$L %*% t(out$L), Sigma, tolerance = 1e-10)
    expect_equal(out$L[1, 2], 0, tolerance = 1e-12)
  })
})

describe(".dd_choice_descriptive_eta", {
  it("eta = Z theta + Zre (L b) per subject; magnitude scale-invariance holds", {
    # rows 1 & 2 share the SAME ll/ss ratio (65/40 == 130/80) and SAME delay =>
    # identical design rows => identical eta (scale-invariance). Row 3 has a
    # different delay so the design is full rank (the strict rank guard passes).
    dat <- data.frame(
      id = c("a", "a", "a"),
      ss_amount = c(40, 80, 40), ll_amount = c(65, 130, 65),
      delay = c(25, 25, 60), choice = c(1, 1, 0)
    )
    d <- .dd_choice_descriptive_design(dat, predictors = NULL, random_slopes = TRUE)
    theta <- c(1.2, -0.3)
    b <- matrix(c(0.1, -0.2), nrow = 1)
    L <- .dd_chol_sigma(log(c(0.5, 0.5)), 0)$L
    eta <- .dd_choice_descriptive_eta(d$Z, theta, d$Zre, b, L,
                                      subject_id = c(0L, 0L, 0L))
    expect_equal(eta[1], eta[2], tolerance = 1e-10)   # scale-invariant pair
  })
})
