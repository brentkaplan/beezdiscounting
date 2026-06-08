# Compile gate: the C++ ChoiceDiscounting -nll (structural mode) must equal the
# R binomial nll (dbinom, logit link) + the N(0,1) RE prior, to 1e-8, with the
# random intercept zeroed (u = 0). dbinom_robust(y, 1, logit_p) == standard
# binomial at p = plogis(eta), so the R-side reference is plain dbinom.

.choice_gate_setup <- function(equation, beta0, log_gamma, has_intercept,
                               beta_k = log(0.02)) {
  ss <- c(40, 55, 31, 14, 47, 25)
  ll <- c(65, 75, 85, 25, 60, 30)
  delay <- c(25, 61, 14, 19, 160, 7)
  y <- c(0, 1, 1, 0, 1, 0)
  n_obs <- length(y)
  eqn_type <- if (equation == "mazur") 0L else 1L
  k <- exp(beta_k)
  D <- if (equation == "mazur") 1 / (1 + k * delay) else exp(-k * delay)
  gamma <- exp(log_gamma)
  eta <- (if (has_intercept) beta0 else 0) + gamma * ((ll / ss) * D - 1)
  list(
    data = list(
      model = "ChoiceDiscounting", mode = 0L, eqn_type = eqn_type,
      has_intercept = as.integer(has_intercept),
      choice = y, ss_amount = ss, ll_amount = ll, delay = delay,
      subject_id = rep(0L, n_obs),
      X = matrix(1, nrow = n_obs, ncol = 1L),
      n_obs = n_obs, n_subjects = 1L
    ),
    parameters = list(
      beta_k = beta_k, log_sigma_u = log(0.5), log_gamma = log_gamma,
      beta0 = beta0, u = matrix(0, nrow = 1L, ncol = 1L)
    ),
    eta = eta, y = y
  )
}

.choice_gate_obj <- function(g, has_intercept) {
  map <- if (isTRUE(has_intercept)) NULL else list(beta0 = factor(NA))
  TMB::MakeADFun(data = g$data, parameters = g$parameters, map = map,
                 DLL = "beezdiscounting", silent = TRUE)
}

.choice_gate_nll <- function(g, has_intercept) {
  obj <- .choice_gate_obj(g, has_intercept)
  as.numeric(obj$fn(obj$par))
}

describe("ChoiceDiscounting compile gate (structural)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("links the ChoiceDiscounting template via the package DLL", {
    expect_true("beezdiscounting" %in% names(getLoadedDLLs()))
  })

  it("C++ -nll == R binomial nll to 1e-8 (mazur, intercept off)", {
    g <- .choice_gate_setup("mazur", beta0 = 0, log_gamma = log(3),
                            has_intercept = FALSE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    obj <- .choice_gate_obj(g, FALSE)
    expect_false("beta0" %in% names(obj$par))
    expect_equal(as.numeric(obj$fn(obj$par)), r_nll, tolerance = 1e-8)
  })

  it("C++ -nll == R binomial nll to 1e-8 (mazur, intercept on)", {
    g <- .choice_gate_setup("mazur", beta0 = 0.4, log_gamma = log(3),
                            has_intercept = TRUE)
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    obj <- .choice_gate_obj(g, TRUE)
    expect_true("beta0" %in% names(obj$par))
    expect_equal(as.numeric(obj$fn(obj$par)), r_nll, tolerance = 1e-8)
  })

  it("C++ -nll == R binomial nll to 1e-8 (exponential, intercept off)", {
    g <- .choice_gate_setup("exponential", beta0 = 0, log_gamma = log(2),
                            has_intercept = FALSE, beta_k = log(5e-3))
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    expect_equal(.choice_gate_nll(g, FALSE), r_nll, tolerance = 1e-8)
  })

  it("matches R dbinom under saturated logits (|eta| > 25) and stays finite", {
    # large log_gamma drives some rows' eta past +/-25 (saturated logits)
    g <- .choice_gate_setup("mazur", beta0 = 0, log_gamma = log(28),
                            has_intercept = FALSE, beta_k = log(0.02))
    expect_gt(max(abs(g$eta)), 25)         # confirm we truly exercised saturation
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    cpp_nll <- .choice_gate_nll(g, FALSE)
    expect_true(is.finite(cpp_nll))
    expect_equal(cpp_nll, r_nll, tolerance = 1e-8)
  })

  it("mode = 1 (descriptive) errors — not implemented in this slice", {
    g <- .choice_gate_setup("mazur", 0, log(3), FALSE)
    g$data$mode <- 1L
    expect_error(
      .choice_gate_obj(g, FALSE),
      regexp = "descriptive|mode|not implemented|Unknown"
    )
  })
})
