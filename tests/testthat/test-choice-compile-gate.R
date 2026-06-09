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
      Z = matrix(0, n_obs, 1L), Zre = matrix(0, n_obs, 1L), n_re = 0L,
      n_obs = n_obs, n_subjects = 1L
    ),
    parameters = list(
      beta_k = beta_k, log_sigma_u = log(0.5), log_gamma = log_gamma,
      beta0 = beta0, u = matrix(0, nrow = 1L, ncol = 1L),
      theta = 0, log_sd_re = rep(log(0.5), 2L), cor_re = 0,
      b = matrix(0, nrow = 1L, ncol = 2L)
    ),
    eta = eta, y = y
  )
}

.choice_gate_obj <- function(g, has_intercept) {
  map <- list(
    theta = factor(rep(NA, length(g$parameters$theta))),
    log_sd_re = factor(rep(NA, length(g$parameters$log_sd_re))),
    cor_re = factor(rep(NA, length(g$parameters$cor_re))),
    b = factor(rep(NA, length(g$parameters$b)))
  )
  if (!isTRUE(has_intercept)) map$beta0 <- factor(NA)
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

  it("mode = 2 (unknown) errors — only 0/1 are implemented", {
    g <- .choice_gate_setup("mazur", 0, log(3), FALSE)
    g$data$mode <- 2L
    expect_error(
      .choice_gate_obj(g, FALSE),
      regexp = "mode|Unknown|implemented"
    )
  })
})

# --- Descriptive (mode 1) compile gate ---------------------------------------
.choice_desc_gate_setup <- function(theta, log_sd_re, cor_re, b_vals,
                                     random_slopes = TRUE) {
  dat <- data.frame(
    id = rep(c("1", "2"), each = 3),
    ss_amount = c(40, 55, 31, 14, 47, 25),
    ll_amount = c(65, 75, 85, 25, 60, 30),
    delay = c(25, 61, 14, 19, 160, 7),
    choice = c(0, 1, 1, 0, 1, 0)
  )
  d <- .dd_choice_descriptive_design(dat, predictors = NULL,
                                     random_slopes = random_slopes)
  n_obs <- nrow(dat)
  n_subj <- 2L
  subject_id <- rep(0:1, each = 3)
  q <- d$q
  b <- matrix(b_vals, nrow = n_subj, ncol = max(q, 1L))
  L <- if (q > 0L) .dd_chol_sigma(log_sd_re, cor_re)$L else matrix(0, 1, 1)
  eta <- .dd_choice_descriptive_eta(
    d$Z, theta, d$Zre, b[, seq_len(max(q, 0L)), drop = FALSE], L, subject_id)
  list(
    data = list(
      model = "ChoiceDiscounting", mode = 1L, eqn_type = 0L, has_intercept = 0L,
      choice = dat$choice, ss_amount = dat$ss_amount, ll_amount = dat$ll_amount,
      delay = dat$delay, subject_id = subject_id,
      X = matrix(1, n_obs, 1L), Z = d$Z, Zre = d$Zre,
      n_obs = n_obs, n_subjects = n_subj, n_re = q
    ),
    parameters = list(
      beta_k = 0, log_sigma_u = log(0.5), log_gamma = 0, beta0 = 0,
      u = matrix(0, n_subj, 1L),
      theta = theta, log_sd_re = log_sd_re, cor_re = cor_re,
      b = matrix(b_vals, nrow = n_subj, ncol = max(q, 1L))
    ),
    eta = eta, y = dat$choice, q = q, b = b
  )
}

.choice_desc_gate_obj <- function(g) {
  map <- list(beta_k = factor(NA), log_sigma_u = factor(NA),
              log_gamma = factor(NA), beta0 = factor(NA),
              u = factor(rep(NA, length(g$parameters$u))))
  if (g$q == 0L) {
    map$log_sd_re <- factor(rep(NA, length(g$parameters$log_sd_re)))
    map$cor_re <- factor(rep(NA, length(g$parameters$cor_re)))
    map$b <- factor(rep(NA, length(g$parameters$b)))
  }
  TMB::MakeADFun(data = g$data, parameters = g$parameters, map = map,
                 DLL = "beezdiscounting", silent = TRUE)
}

describe("ChoiceDiscounting compile gate (descriptive)", {
  skip_on_cran(); skip_if_not_installed("TMB")

  it("C++ -nll == R binomial nll + N(0,1) prior to 1e-8 (q=2, b=0)", {
    g <- .choice_desc_gate_setup(theta = c(1.5, -0.4),
                                 log_sd_re = log(c(0.6, 0.3)), cor_re = atanh(-0.3),
                                 b_vals = c(0, 0, 0, 0))
    re_prior <- -sum(dnorm(g$b, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    obj <- .choice_desc_gate_obj(g)
    expect_equal(as.numeric(obj$fn(obj$par)), r_nll, tolerance = 1e-8)
  })

  it("C++ -nll == R binomial nll + N(0,1) prior to 1e-8 (q=2, b!=0)", {
    g <- .choice_desc_gate_setup(theta = c(1.2, -0.5),
                                 log_sd_re = log(c(0.8, 0.5)), cor_re = atanh(0.4),
                                 b_vals = c(0.7, -0.6, -0.2, 0.9))
    re_prior <- -sum(dnorm(g$b, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    obj <- .choice_desc_gate_obj(g)
    expect_equal(as.numeric(obj$fn(obj$par)), r_nll, tolerance = 1e-8)
  })

  it("pooled (q=0, random_slopes=FALSE): C++ -nll == R binomial nll, no RE prior", {
    g <- .choice_desc_gate_setup(theta = c(1.0, -0.3),
                                 log_sd_re = log(c(0.5, 0.5)), cor_re = 0,
                                 b_vals = 0, random_slopes = FALSE)
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE))  # no RE prior
    obj <- .choice_desc_gate_obj(g)
    expect_equal(as.numeric(obj$fn(obj$par)), r_nll, tolerance = 1e-8)
  })

  it("matches R dbinom under saturated descriptive logits and stays finite", {
    g <- .choice_desc_gate_setup(theta = c(40, -40),
                                 log_sd_re = log(c(0.5, 0.5)), cor_re = 0,
                                 b_vals = c(0, 0, 0, 0))
    expect_gt(max(abs(g$eta)), 25)
    re_prior <- -sum(dnorm(g$b, 0, 1, log = TRUE))
    r_nll <- -sum(dbinom(g$y, 1, plogis(g$eta), log = TRUE)) + re_prior
    cpp <- as.numeric(.choice_desc_gate_obj(g)$fn(.choice_desc_gate_obj(g)$par))
    expect_true(is.finite(cpp))
    expect_equal(cpp, r_nll, tolerance = 1e-8)
  })
})
