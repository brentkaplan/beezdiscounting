describe("simulate_dd_choice() structural", {
  it("is exported and returns the generic choice frame", {
    expect_true("simulate_dd_choice" %in% getNamespaceExports("beezdiscounting"))
    sim <- simulate_dd_choice(n_subjects = 10, seed = 1)
    expect_s3_class(sim, "tbl_df")
    expect_true(all(c("id", "ss_amount", "ll_amount", "delay", "choice") %in% names(sim)))
    expect_true(all(sim$choice %in% c(0, 1)))
  })
  it("P(LL) tracks the structural curve in expectation (sigma_u ~ 0)", {
    set.seed(7)
    sim <- simulate_dd_choice(n_subjects = 400, log_k_pop = log(0.02),
                              sigma_u = 1e-6, gamma = 4, beta0 = 0,
                              equation = "mazur", seed = 7)
    k <- 0.02; gamma <- 4
    D <- 1 / (1 + k * sim$delay)
    p_true <- plogis(gamma * ((sim$ll_amount / sim$ss_amount) * D - 1))
    by_bin <- tapply(sim$choice, cut(p_true, 5), mean)
    centers <- tapply(p_true, cut(p_true, 5), mean)
    expect_lt(max(abs(by_bin - centers), na.rm = TRUE), 0.07)
  })
})

describe("simulate_dd_choice (descriptive)", {
  it("returns the canonical frame and is seed-deterministic", {
    a <- simulate_dd_choice(n_subjects = 15, mode = "descriptive",
                            theta = c(1.5, -0.4), re_sd = c(0.5, 0.3),
                            re_cor = -0.2, seed = 7)
    b <- simulate_dd_choice(n_subjects = 15, mode = "descriptive",
                            theta = c(1.5, -0.4), re_sd = c(0.5, 0.3),
                            re_cor = -0.2, seed = 7)
    expect_named(a, c("id", "ss_amount", "ll_amount", "delay", "choice"))
    expect_true(all(a$choice %in% c(0, 1)))
    expect_identical(a, b)
  })

  it("empirical per-subject slope (co)variance approximates Sigma (loose)", {
    sim <- simulate_dd_choice(n_subjects = 400, mode = "descriptive",
                              theta = c(1.5, -0.4), re_sd = c(0.6, 0.35),
                              re_cor = -0.4, seed = 3, return_truth = TRUE)
    b <- attr(sim, "subject_slopes")          # n_subjects x 2 realized slopes
    emp <- stats::cov(b)
    expect_equal(sqrt(emp[1, 1]), 0.6, tolerance = 0.12)
    expect_equal(sqrt(emp[2, 2]), 0.35, tolerance = 0.12)
    expect_equal(emp[1, 2] / sqrt(emp[1, 1] * emp[2, 2]), -0.4, tolerance = 0.2)
  })
})

describe("simulate_dd_choice() RNG hygiene", {
  it("restores the caller's RNG state when seed is supplied", {
    skip_on_cran()
    set.seed(999)
    before <- get(".Random.seed", envir = globalenv())
    invisible(simulate_dd_choice(n_subjects = 6, seed = 1))
    expect_identical(get(".Random.seed", envir = globalenv()), before)

    set.seed(999)
    expected <- runif(1)
    set.seed(999)
    invisible(simulate_dd_choice(n_subjects = 6, seed = 1))
    expect_identical(runif(1), expected)
  })

  it("simulated output for a given seed is unchanged by the RNG-restore fix", {
    sim <- simulate_dd_choice(n_subjects = 6, seed = 1)
    # Value computed BEFORE the RNG-restore fix (from git HEAD's
    # simulate_dd_choice(), which called set.seed(seed) without saving/
    # restoring the caller's state) -- pinned so the fix cannot alter
    # simulated output.
    expect_identical(sim$choice[1:6], c(0L, 0L, 1L, 1L, 0L, 0L))
  })
})
