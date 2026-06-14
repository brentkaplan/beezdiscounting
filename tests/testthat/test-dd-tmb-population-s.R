# tests/testthat/test-dd-tmb-population-s.R
# Population-level (random-effects-zero) discounting exponent s used by predict()
# at level = "population". For an s-target 2-RE fit the kernel maps log_s through
# the per-subject SOFT clamp (re = 0 for the population curve), so the population s
# must be .dd_soft_clamp_s_log(log_s), NOT exp(log_s) -- they differ near the
# [0.05, 20] bounds. (Regression guard for the Codex capstone 2026-06-14 finding.)

describe(".dd_tmb_population_s", {
  it("returns 1 when the equation has no s (mazur/exponential)", {
    expect_equal(.dd_tmb_population_s(has_s = FALSE, is_s_re = FALSE, log_s = NA_real_), 1)
  })

  it("returns exp(log_s) for a population-only s (1-RE / phi-target GM/Rachlin)", {
    # The kernel uses Type s = exp(log_s) on those paths (no per-subject soft clamp),
    # so predict must match exp(log_s) even near a bound.
    expect_equal(.dd_tmb_population_s(has_s = TRUE, is_s_re = FALSE, log_s = log(1.4)), 1.4)
    expect_equal(.dd_tmb_population_s(has_s = TRUE, is_s_re = FALSE, log_s = log(0.05) - 1),
                 exp(log(0.05) - 1))
  })

  it("returns the SOFT-clamped population s for an s-target 2-RE fit", {
    # Interior: matches both exp and the soft clamp (they coincide).
    expect_equal(.dd_tmb_population_s(has_s = TRUE, is_s_re = TRUE, log_s = log(1.4)),
                 .dd_soft_clamp_s_log(log(1.4)))
    # Near/below the lower bound the soft clamp DIVERGES from exp(log_s); the
    # population predict must use the soft clamp (kernel-consistent), not exp.
    lb <- log(0.05) - 1
    expect_equal(.dd_tmb_population_s(has_s = TRUE, is_s_re = TRUE, log_s = lb),
                 .dd_soft_clamp_s_log(lb))
    expect_gt(abs(.dd_tmb_population_s(has_s = TRUE, is_s_re = TRUE, log_s = lb) - exp(lb)), 1e-3)
  })
})
