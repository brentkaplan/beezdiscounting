# tests/testthat/test-tmb-compile-gate.R
# Compile gate: the C++ MixedDiscounting -nll must equal the verified R SLT
# log-density (and Gaussian density) to 1e-8, isolating the observation
# term by zeroing the random effects.
#
# The R-side density is the canonical .dd_slt_logpdf from R/dd-density.R
# (Task F.1), made available by the `devtools::load_all()` in the run
# commands below. Do NOT inline a local copy here — there is a single source
# of truth for the SLT density. (Gaussian compares against dnorm(log=TRUE).)
#
# *** Marginal vs conditional note ***
# MakeADFun is built WITHOUT `random =` so that obj$fn(par) returns the full
# JOINT (conditional) negative log-likelihood:
#   sum_i  -log f_obs(y_i | mu_i, aux)  +  sum_j  -dnorm(u_j, 0, 1, log=TRUE)
# Using random = "u" would return the Laplace-marginal nll (u integrated out),
# which does NOT equal the R conditional density sum even at u = 0.  The
# RE prior is accounted for explicitly on the R side (re_prior below).

# Build a single-subject, intercept-only TMB data/parameter set and the
# matching R-side prediction for a given equation/family.
.gate_setup <- function(eqn_type, family, beta0, aux) {
  delays <- c(7, 30, 180, 365, 730, 1460, 2920)
  y <- c(0.95, 0.80, 0.50, 0.35, 0.20, 0.08, 0.00)
  n_obs <- length(delays)
  k <- exp(beta0)
  mu_raw <- if (eqn_type == 0L) 1 / (1 + k * delays) else exp(-k * delays)
  mu <- pmin(pmax(mu_raw, 1e-6), 1 - 1e-6)
  list(
    data = list(
      model = "MixedDiscounting",
      y = y, x = delays,
      subject_id = rep(0L, n_obs),
      X = matrix(1, nrow = n_obs, ncol = 1L),
      eqn_type = eqn_type, family = family,
      n_obs = n_obs, n_subjects = 1L
    ),
    parameters = list(
      beta_k = beta0,
      log_sigma_u = log(0.5),     # arbitrary; RE drops out because u = 0
      log_aux = log(aux),
      u = matrix(0, nrow = 1L, ncol = 1L)
    ),
    mu = mu, y = y
  )
}

# Flatten parameters list to a named numeric vector matching TMB's par order.
.gate_par <- function(g) {
  c(beta_k = g$parameters$beta_k,
    log_sigma_u = g$parameters$log_sigma_u,
    log_aux = g$parameters$log_aux,
    u = 0)
}

describe("MixedDiscounting compile gate", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  it("compiles and links the MixedDiscounting template", {
    # Loading the installed package DLL is enough; if devtools::load_all was
    # used, the DLL is already loaded. getLoadedDLLs() must list beezdiscounting.
    expect_true("beezdiscounting" %in% names(getLoadedDLLs()))
  })

  it("C++ -nll matches the R SLT log-density to 1e-8 (mazur, sltb)", {
    g <- .gate_setup(eqn_type = 0L, family = 0L, beta0 = log(0.01), aux = 8)
    # No random = "u": obj$fn returns the JOINT (conditional) nll so the
    # comparison against the R obs density + RE prior is exact.
    obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                          DLL = "beezdiscounting", silent = TRUE)
    cpp_nll <- obj$fn(.gate_par(g))
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))           # one subject, u = 0
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches the R SLT log-density to 1e-8 (exponential, sltb)", {
    g <- .gate_setup(eqn_type = 1L, family = 0L, beta0 = log(2e-4), aux = 8)
    obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                          DLL = "beezdiscounting", silent = TRUE)
    cpp_nll <- obj$fn(.gate_par(g))
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-.dd_slt_logpdf(g$y, g$mu, 8)) + re_prior
    expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
  })

  it("C++ -nll matches the R Gaussian density to 1e-8 (mazur, gaussian)", {
    g <- .gate_setup(eqn_type = 0L, family = 1L, beta0 = log(0.01), aux = 0.1)
    obj <- TMB::MakeADFun(data = g$data, parameters = g$parameters,
                          DLL = "beezdiscounting", silent = TRUE)
    cpp_nll <- obj$fn(.gate_par(g))
    re_prior <- -sum(dnorm(0, 0, 1, log = TRUE))
    r_nll <- sum(-dnorm(g$y, g$mu, 0.1, log = TRUE)) + re_prior
    expect_equal(as.numeric(cpp_nll), r_nll, tolerance = 1e-8)
  })
})
