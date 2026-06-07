# B3: the fixed-effect design must reject rank-deficient (empty-cell / aliased)
# interaction designs rather than silently fitting a non-estimable model whose
# EMMs would average over a phantom cell.

describe("design rank guard (B3)", {
  skip_on_cran()
  skip_if_not_installed("TMB")

  # Build a 2x2 between-subject design. When omit_a2b2 = TRUE the (a2, b2) cell
  # has no subjects, so model.matrix(~ A * B) has an all-zero Aa2:Bb2 column.
  .cell_data <- function(omit_a2b2 = TRUE, seed = 3) {
    sim <- simulate_dd_ip(n_subjects = 20, family = "gaussian", seed = seed)
    ids <- levels(sim$id)
    if (omit_a2b2) {
      a <- rep(c("a1", "a1", "a2"), length.out = length(ids))   # 3 cells only
      b <- rep(c("b1", "b2", "b1"), length.out = length(ids))
    } else {
      a <- rep(c("a1", "a1", "a2", "a2"), length.out = length(ids))  # all 4 cells
      b <- rep(c("b1", "b2", "b1", "b2"), length.out = length(ids))
    }
    A_by <- stats::setNames(a, ids)
    B_by <- stats::setNames(b, ids)
    sim$A <- factor(A_by[as.character(sim$id)])
    sim$B <- factor(B_by[as.character(sim$id)])
    sim
  }

  it("aborts on a rank-deficient interaction design (empty a2:b2 cell)", {
    dat <- .cell_data(omit_a2b2 = TRUE)
    expect_error(
      fit_dd_tmb(dat, family = "gaussian", factors = c("A", "B"),
                 factor_interaction = TRUE, multi_start = FALSE, verbose = 0),
      "rank-deficient|non-estimable|aliased"
    )
  })

  it("the same empty-cell data fits under an ADDITIVE design (main effects identified)", {
    dat <- .cell_data(omit_a2b2 = TRUE)
    expect_no_error(
      fit_dd_tmb(dat, family = "gaussian", factors = c("A", "B"),
                 factor_interaction = FALSE, multi_start = FALSE, verbose = 0)
    )
  })

  it("a full 2x2 interaction design fits (no false positive)", {
    dat <- .cell_data(omit_a2b2 = FALSE)
    expect_no_error(
      fit_dd_tmb(dat, family = "gaussian", factors = c("A", "B"),
                 factor_interaction = TRUE, multi_start = FALSE, verbose = 0)
    )
  })
})
