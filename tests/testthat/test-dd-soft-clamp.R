# tests/testthat/test-dd-soft-clamp.R
# Pure-R unit tests for the two-sided softplus soft clamp used (mirrored) by the
# s-target 2-RE kernel. C++ parity is covered separately by the compile gate.

describe(".dd_soft_clamp_s_log", {
  it("is bounded inside (0.05, 20), monotone, and near-identity in the interior", {
    u <- seq(-12, 12, length.out = 2001)          # s_raw = exp(u) far past both bounds
    s <- .dd_soft_clamp_s_log(u)
    expect_true(all(s >= 0.05 - 1e-9 & s <= 20 + 1e-9))  # closed bound [0.05, 20] up to float
                                                   # roundoff at saturation (exp() rounds the
                                                   # bound by ~1 ULP for |u| >> the bounds;
                                                   # strictly inside for any finite fit latent)
    expect_true(all(diff(s) >= -1e-12))            # monotone up to float roundoff in the
                                                   # saturated tail (Codex finding #1)
    ui <- log(c(0.2, 0.5, 1, 1.4, 3, 5))           # deep interior
    expect_lt(max(abs(.dd_soft_clamp_s_log(ui) / exp(ui) - 1)), 1e-3)  # ~identity
  })

  it("has a strictly positive gradient where the HARD clamp was flat (=> no hang)", {
    # Evaluate INSIDE the old flat zones (0.1 past each bound) so the central
    # difference does not straddle the interior side. There the hard clamp's
    # gradient is exactly 0; the soft clamp's must be > 0 (Codex finding #2).
    h <- 1e-5
    for (bd in c(.dd_s_log_lower - 0.1, .dd_s_log_upper + 0.1)) {
      g <- (.dd_soft_clamp_s_log(bd + h) - .dd_soft_clamp_s_log(bd - h)) / (2 * h)
      expect_gt(g, 1e-5)
    }
  })

  it("differs from the hard clamp at the bounds, agrees in the deep interior", {
    hard <- function(u) pmin(pmax(exp(u), 0.05), 20)
    for (u in c(.dd_s_log_lower, .dd_s_log_upper)) {           # binding points: differ
      expect_gt(abs(.dd_soft_clamp_s_log(u) - hard(u)), 1e-4)
    }
    for (u in log(c(0.5, 1, 1.4, 5))) {                        # deep interior: identical
      expect_lt(abs(.dd_soft_clamp_s_log(u) - hard(u)), 1e-8)
    }
  })
})
