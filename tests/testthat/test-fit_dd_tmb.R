describe(".dd_tmb_prepare_data()", {
  it("0-indexes subject_id and aligns vectors with subject_levels", {
    dat <- data.frame(
      id = factor(c("b", "b", "a", "a"), levels = c("b", "a")),
      x = c(7, 30, 7, 30),
      y = c(0.9, 0.5, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$subject_levels, c("a", "b"))
    expect_equal(prep$n_subjects, 2L)
    expect_equal(prep$n_obs, 4L)
    expect_equal(min(prep$subject_id), 0L)
    expect_equal(max(prep$subject_id), 1L)
    # subject_id maps each row to the sorted-level index
    expect_equal(prep$subject_id, c(1L, 1L, 0L, 0L))
    expect_type(prep$subject_id, "integer")
  })

  it("drops rows with NA in y/x/id and keeps x,y,subject_id row-coherent", {
    dat <- data.frame(
      id = c("a", "a", "b", "b"),
      x = c(7, 30, 7, NA),
      y = c(0.9, NA, 0.8, 0.4)
    )
    prep <- .dd_tmb_prepare_data(dat, y_var = "y", x_var = "x", id_var = "id")
    expect_equal(prep$n_obs, 2L)
    expect_equal(prep$y, c(0.9, 0.8))
    expect_equal(prep$x, c(7, 7))
    expect_equal(nrow(prep$data), 2L)
    # surviving rows keep their parallel arrays aligned
    expect_equal(length(prep$x), length(prep$subject_id))
  })

  it("renames caller columns to canonical id/x/y in $data", {
    dat <- data.frame(pid = "p1", delay = 7, ip = 0.5)
    prep <- .dd_tmb_prepare_data(dat, y_var = "ip", x_var = "delay", id_var = "pid")
    expect_true(all(c("id", "x", "y") %in% names(prep$data)))
    expect_equal(names(prep$data)[1:3], c("id", "x", "y"))
  })

  it("retains extra cols and complete-cases id/x/y + extras ONCE (row-coherent)", {
    dat <- data.frame(
      id = c("a", "a", "b", "b"),
      x = c(7, 30, 7, 30),
      y = c(0.9, 0.5, 0.8, 0.4),
      grp = factor(c("ctrl", "ctrl", NA, "trt")),  # one NA in an extra col
      age = c(20, 20, 30, 30),
      stringsAsFactors = FALSE
    )
    prep <- .dd_tmb_prepare_data(dat, "y", "x", "id",
      extra_cols = c("grp", "age"))
    # the row with NA grp is dropped because grp is a modeling column
    expect_equal(prep$n_obs, 3L)
    # $data keeps id/x/y first then the retained extras, all row-aligned
    expect_equal(names(prep$data), c("id", "x", "y", "grp", "age"))
    expect_equal(prep$data$y, c(0.9, 0.5, 0.4))
    expect_equal(as.character(prep$data$grp), c("ctrl", "ctrl", "trt"))
    # parallel arrays match $data exactly (one filtered row order)
    expect_equal(prep$y, prep$data$y)
    expect_equal(prep$x, prep$data$x)
    expect_equal(length(prep$subject_id), nrow(prep$data))
    # unused factor levels dropped on the surviving frame
    expect_setequal(levels(prep$data$grp), c("ctrl", "trt"))
  })
})

describe(".dd_tmb_build_design()", {
  it("returns an intercept-only X (+ rhs formula + contrasts) with no factors", {
    dat <- data.frame(id = "a", x = 7, y = 0.5)
    d <- .dd_tmb_build_design(dat)
    expect_equal(colnames(d$X), "(Intercept)")
    expect_equal(ncol(d$X), 1L)
    # rhs is a one-sided formula from build_fixed_rhs(), not a string
    expect_s3_class(d$rhs, "formula")
    expect_equal(rlang::f_text(d$rhs), "1")
    # contrasts attr is carried through (NULL when no factors)
    expect_true("contrasts" %in% names(d))
  })

  it("expands a single between-subject factor into contrast columns", {
    dat <- data.frame(
      id = rep(c("a", "b", "c", "d"), each = 1),
      x = 7,
      y = 0.5,
      grp = factor(c("ctrl", "ctrl", "trt", "trt"))
    )
    d <- .dd_tmb_build_design(dat, factors = "grp")
    expect_equal(ncol(d$X), 2L)
    expect_true("grptrt" %in% colnames(d$X))
    expect_match(rlang::f_text(d$rhs), "grp")
    # the per-factor contrast scheme is stored for the emmeans ref grid (E.1)
    expect_true("grp" %in% names(d$contrasts))
  })

  it("adds an interaction term when factor_interaction = TRUE", {
    dat <- data.frame(
      id = letters[1:4], x = 7, y = 0.5,
      a = factor(c("x", "x", "z", "z")),
      b = factor(c("p", "q", "p", "q"))
    )
    d <- .dd_tmb_build_design(dat, factors = c("a", "b"),
                              factor_interaction = TRUE)
    expect_true(any(grepl(":", colnames(d$X))))
  })

  it("adds continuous covariate main-effect columns", {
    dat <- data.frame(id = letters[1:3], x = 7, y = 0.5, age = c(20, 30, 40))
    d <- .dd_tmb_build_design(dat, continuous_covariates = "age")
    expect_true("age" %in% colnames(d$X))
  })
})

describe(".dd_tmb_build_tmb_data()", {
  it("assembles the TMB data list with the contract names and enum ints", {
    prep <- list(
      y = c(0.9, 0.5), x = c(7, 30), subject_id = c(0L, 0L),
      subject_levels = "a", n_subjects = 1L, n_obs = 2L,
      data = data.frame(id = "a", x = c(7, 30), y = c(0.9, 0.5))
    )
    design <- list(X = matrix(1, nrow = 2, ncol = 1,
                              dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "mazur",
                                 family = "sltb")
    expect_equal(td$model, "MixedDiscounting")
    expect_equal(td$eqn_type, 0L)   # mazur
    expect_equal(td$family, 0L)     # sltb
    expect_equal(td$n_obs, 2L)
    expect_equal(td$n_subjects, 1L)
    expect_equal(td$subject_id, c(0L, 0L))
    expect_true(is.matrix(td$X))
  })

  it("maps exponential -> 1 and gaussian -> 1", {
    prep <- list(y = 0.5, x = 7, subject_id = 0L, subject_levels = "a",
                 n_subjects = 1L, n_obs = 1L,
                 data = data.frame(id = "a", x = 7, y = 0.5))
    design <- list(X = matrix(1, 1, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    td <- .dd_tmb_build_tmb_data(prep, design, equation = "exponential",
                                 family = "gaussian")
    expect_equal(td$eqn_type, 1L)
    expect_equal(td$family, 1L)
  })
})

describe(".dd_tmb_default_starts()", {
  it("derives log_k intercept from median y at min delay (mazur)", {
    # median y at min delay = 0.8 -> mu=0.8 -> k = (1/mu - 1)/x_min
    prep <- list(
      y = c(0.8, 0.8, 0.3, 0.3), x = c(7, 7, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(7, 7, 365, 365), y = c(0.8, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 1L)
    k_implied <- (1 / 0.8 - 1) / 7
    expect_equal(st$beta_k[1], log(k_implied), tolerance = 1e-8)
    expect_equal(st$log_sigma_u, log(0.5))
    expect_equal(st$log_aux, log(8))          # sltb
    expect_equal(dim(st$u), c(2L, 1L))
  })

  it("uses log(0.1) for log_aux under gaussian", {
    prep <- list(y = c(0.8, 0.3), x = c(7, 365), subject_id = c(0L, 0L),
                 subject_levels = "a", n_subjects = 1L, n_obs = 2L,
                 data = data.frame(id = "a", x = c(7, 365), y = c(0.8, 0.3)))
    design <- list(X = matrix(1, 2, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "gaussian")
    expect_equal(st$log_aux, log(0.1))
  })

  it("derives log_k intercept from exponential inversion (equation = 'exponential')", {
    # median y at min delay (x=7) = 0.8 -> mu=0.8 -> k = -log(mu)/x_min
    prep <- list(
      y = c(0.8, 0.8, 0.3, 0.3), x = c(7, 7, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(7, 7, 365, 365), y = c(0.8, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                 equation = "exponential")
    k_implied <- -log(0.8) / 7
    expect_equal(st$beta_k[1], log(k_implied), tolerance = 1e-8)
  })

  it("falls back to log(0.01) intercept when x_min <= 0", {
    # x_min = 0: inversion would divide by zero; must use safe default k0 = 0.01
    prep <- list(
      y = c(0.9, 0.8, 0.3, 0.3), x = c(0, 0, 365, 365),
      subject_id = c(0L, 1L, 0L, 1L), subject_levels = c("a", "b"),
      n_subjects = 2L, n_obs = 4L,
      data = data.frame(id = c("a", "b", "a", "b"),
                        x = c(0, 0, 365, 365), y = c(0.9, 0.8, 0.3, 0.3))
    )
    design <- list(X = matrix(1, 4, 1, dimnames = list(NULL, "(Intercept)")),
                   rhs = "~ 1")
    st_mazur <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                       equation = "mazur")
    st_exp <- .dd_tmb_default_starts(prep, design, family = "sltb",
                                     equation = "exponential")
    expect_equal(st_mazur$beta_k[1], log(0.01), tolerance = 1e-8)
    expect_equal(st_exp$beta_k[1],   log(0.01), tolerance = 1e-8)
  })

  it("zero-pads beta_k for multi-column designs", {
    prep <- list(y = c(0.8, 0.5), x = c(7, 7), subject_id = c(0L, 1L),
                 subject_levels = c("a", "b"), n_subjects = 2L, n_obs = 2L,
                 data = data.frame(id = c("a", "b"), x = 7, y = c(0.8, 0.5)))
    design <- list(
      X = matrix(c(1, 1, 0, 1), 2, 2,
                 dimnames = list(NULL, c("(Intercept)", "grptrt"))),
      rhs = "~ grp")
    st <- .dd_tmb_default_starts(prep, design, family = "sltb")
    expect_length(st$beta_k, 2L)
    expect_equal(st$beta_k[2], 0)
  })
})
