describe(".dd_choice_structural_eta()", {
  it("computes beta0 + gamma*((ll/ss)*D(k,delay) - 1) for mazur", {
    k <- 0.02; gamma <- 3; beta0 <- 0.1
    ss <- c(40, 31); ll <- c(65, 85); delay <- c(25, 14)
    D <- 1 / (1 + k * delay)                       # mazur
    expected <- beta0 + gamma * ((ll / ss) * D - 1)
    got <- .dd_choice_structural_eta(k = k, ss_amount = ss, ll_amount = ll,
                                     delay = delay, equation = "mazur",
                                     gamma = gamma, beta0 = beta0)
    expect_equal(got, expected, tolerance = 1e-12)
  })

  it("supports exponential and a zero intercept", {
    k <- 0.01; gamma <- 2
    ss <- 40; ll <- 80; delay <- 30
    D <- exp(-k * delay)
    expected <- 0 + gamma * ((ll / ss) * D - 1)
    got <- .dd_choice_structural_eta(k, ss, ll, delay, "exponential", gamma, 0)
    expect_equal(got, expected, tolerance = 1e-12)
  })

  it("is magnitude-scale-invariant (same ll/ss & delay -> same eta)", {
    k <- 0.02; gamma <- 3; beta0 <- 0
    # two trials, identical ll/ss ratio (2x) and delay, different absolute size
    e_small <- .dd_choice_structural_eta(k, 10, 20, 30, "mazur", gamma, beta0)
    e_big   <- .dd_choice_structural_eta(k, 50, 100, 30, "mazur", gamma, beta0)
    expect_equal(e_small, e_big, tolerance = 1e-12)
  })

  it("beta0 = 0 gives eta = 0 (P = 0.5) exactly at the indifference k", {
    # mazur indifference: ll/ss = 1 + k*delay  =>  (ll/ss)*D = 1  =>  eta = 0
    delay <- 25; k <- 0.02; ss <- 40; ll <- ss * (1 + k * delay)
    eta <- .dd_choice_structural_eta(k, ss, ll, delay, "mazur", gamma = 5, beta0 = 0)
    expect_equal(eta, 0, tolerance = 1e-10)
    expect_equal(stats::plogis(eta), 0.5, tolerance = 1e-10)
  })

  it("errors on an unsupported equation", {
    expect_error(.dd_choice_structural_eta(0.01, 40, 80, 30, "rachlin", 2, 0),
                 "mazur|exponential|unsupported")
  })
})
