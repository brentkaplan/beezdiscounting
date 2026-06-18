# ==============================================================================
# plot() methods for the modern discounting tiers (TMB, choice, brms) and the
# group-comparison object, plus a plot_qq() random-effect normality diagnostic.
#
# All public methods live here (rather than the per-tier *-methods.R files) so
# the plotting surface stays in one place and edits do not churn the larger
# method files. Each method returns a ggplot object; nothing is printed.
#
# Shared internals build a NORMALIZED curve frame -- x, .value, .group, and
# (Bayesian only) .lower/.upper -- so a single assembler can draw every tier
# regardless of whether the underlying predict() returns predict.fixed, .prob,
# or .fitted/.lower/.upper.
# ==============================================================================

#' @importFrom ggplot2 .data
NULL

# Palette: pkgdown primary (population curve) + secondary (bands / accents).
.dd_col_pop <- "#2B4560"
.dd_col_accent <- "#A25F5F"

# Theme shared with plot_dd() for a consistent look across the package.
.dd_plot_theme <- function() {
  beezdemand::theme_apa()
}

# Delay-column name. IP fits store it in param_info$x_var ("x"); the choice
# tiers use the canonical "delay" (and the brms choice fit does not populate
# param_info$x_var), so fall back by class.
.dd_x_var <- function(fit) {
  xv <- fit$param_info$x_var
  if (is.null(xv) || !nzchar(xv)) {
    xv <- if (
      inherits(fit, c("beezdiscounting_choice", "beezdiscounting_choice_brms"))
    ) {
      "delay"
    } else {
      "x"
    }
  }
  xv
}

# Delay-axis sequence over the observed range. On a log10 axis only positive
# delays are representable, so non-positive delays (e.g. the immediate option in
# an IP task) are excluded from the grid.
.dd_x_seq <- function(xvals, n_points, x_trans) {
  xvals <- xvals[is.finite(xvals)]
  if (x_trans == "log10") {
    xvals <- xvals[xvals > 0]
  }
  if (length(unique(xvals)) < 2L) {
    cli::cli_abort(c(
      "Not enough distinct delays to draw a curve.",
      "i" = if (x_trans == "log10") {
        "A {.code log10} delay axis drops non-positive delays; try {.code x_trans = \"linear\"}."
      } else {
        "At least two finite delays are required."
      }
    ))
  }
  r <- range(xvals)
  if (x_trans == "log10") {
    exp(seq(log(r[1]), log(r[2]), length.out = n_points))
  } else {
    seq(r[1], r[2], length.out = n_points)
  }
}

.dd_x_scale <- function(x_trans) {
  if (x_trans == "log10") {
    ggplot2::scale_x_log10()
  } else {
    ggplot2::scale_x_continuous()
  }
}

# Drop non-positive x from a frame on a log10 axis, warning once. Used for the
# observed-point overlay (the curve grid is already positive-only).
.dd_drop_nonpos_x <- function(df, xcol, x_trans, what = "observations") {
  if (x_trans != "log10" || !nrow(df)) {
    return(df)
  }
  keep <- df[[xcol]] > 0 & is.finite(df[[xcol]])
  if (any(!keep)) {
    cli::cli_inform(c(
      "i" = "Dropped {sum(!keep)} {what} at non-positive delay from the {.code log10} axis.",
      " " = "Use {.code x_trans = \"linear\"} to keep them."
    ))
  }
  df[keep, , drop = FALSE]
}

# ------------------------------------------------------------------------------
# Curve-grid construction (factor/covariate aware, via the EMM reference grid).
# ------------------------------------------------------------------------------

# Population newdata: the delay grid crossed with each fitted factor-level combo
# (continuous covariates held at their training mean, or the supplied `at`
# value). Returns the newdata frame plus a `.group` label column. For an
# intercept-only fit there is a single unlabelled population curve.
.dd_curve_newdata <- function(fit, at, n_points, x_trans) {
  x_var <- .dd_x_var(fit)
  xs <- .dd_x_seq(fit$data[[x_var]], n_points, x_trans)

  grid <- .dd_build_emm_ref_grid(fit, at = at)
  if (isTRUE(grid$is_intercept_only)) {
    nd <- data.frame(xs)
    names(nd) <- x_var
    nd$.group <- factor("Population")
    return(nd)
  }

  lc <- grid$level_combos
  nrep <- nrow(lc)
  nd <- lc[rep(seq_len(nrep), each = length(xs)), , drop = FALSE]
  nd[[x_var]] <- rep(xs, times = nrep)

  fac_cols <- intersect(names(lc), grid$use_factors)
  if (length(fac_cols)) {
    lab <- do.call(
      paste,
      c(lapply(fac_cols, function(f) as.character(lc[[f]])), list(sep = ", "))
    )
    nd$.group <- factor(rep(lab, each = length(xs)), levels = unique(lab))
  } else {
    nd$.group <- factor("Population")
  }
  rownames(nd) <- NULL
  nd
}

# Per-subject newdata: each subject's own factor/covariate values (taken from
# the first training row for that id) crossed with the delay grid.
.dd_individual_newdata <- function(fit, ids, n_points, x_trans) {
  x_var <- .dd_x_var(fit)
  id_var <- fit$param_info$id_var
  xs <- .dd_x_seq(fit$data[[x_var]], n_points, x_trans)

  all_ids <- as.character(fit$param_info$subject_levels)
  if (!is.null(ids)) {
    ids <- unique(as.character(ids))
    if (!length(ids)) {
      cli::cli_abort("{.arg ids} is empty.")
    }
    unknown <- setdiff(ids, all_ids)
    if (length(unknown)) {
      cli::cli_abort(c(
        "Unknown {.arg ids}: {.val {unknown}}.",
        "i" = "Fitted subjects: {.val {all_ids}}."
      ))
    }
    all_ids <- ids
  }

  carry <- intersect(
    c(fit$param_info$factors, fit$param_info$continuous_covariates),
    names(fit$data)
  )
  carry <- carry[nzchar(carry) & !is.na(carry)]
  dat <- fit$data
  first_row <- dat[match(all_ids, as.character(dat[[id_var]])), , drop = FALSE]

  nd <- first_row[
    rep(seq_along(all_ids), each = length(xs)),
    c(id_var, carry),
    drop = FALSE
  ]
  nd[[x_var]] <- rep(xs, times = length(all_ids))
  nd$.group <- factor(rep(all_ids, each = length(xs)), levels = all_ids)
  rownames(nd) <- NULL
  nd
}

# ------------------------------------------------------------------------------
# Tier-specific curve values -> normalized frame (x, .value, .group[, .lower,
# .upper]). The IP tiers (tmb, brms) predict the indifference proportion; the
# brms tier additionally carries a credible band.
# ------------------------------------------------------------------------------
.dd_ip_curve_values <- function(fit, nd, level) {
  x_var <- .dd_x_var(fit)
  grp <- nd$.group
  pr <- stats::predict(
    fit,
    newdata = nd[setdiff(names(nd), ".group")],
    type = "response",
    level = if (level == "population") "population" else "subject"
  )
  if (inherits(fit, "beezdiscounting_brms")) {
    out <- data.frame(x = nd[[x_var]], .value = pr$.fitted, .group = grp)
    if (level == "population") {
      out$.lower <- pr$.lower
      out$.upper <- pr$.upper
    }
    return(out)
  }
  value <- if (level == "population") pr$predict.fixed else pr$.fitted
  data.frame(x = nd[[x_var]], .value = value, .group = grp)
}

# ------------------------------------------------------------------------------
# Shared assembler: observed points + faint per-subject lines + population
# line/ribbon. `multi_group` colours one population line per factor combo.
# ------------------------------------------------------------------------------
.dd_assemble_curve <- function(
  curve_pop,
  curve_ind = NULL,
  observed = NULL,
  x_trans = "log10",
  xlab = "Delay",
  ylab = "Indifference point",
  title = NULL
) {
  multi_group <- nlevels(curve_pop$.group) > 1L
  p <- ggplot2::ggplot()

  if (!is.null(observed) && nrow(observed)) {
    p <- p +
      ggplot2::geom_point(
        data = observed,
        ggplot2::aes(x = .data$x, y = .data$y),
        alpha = 0.3,
        size = 1.5,
        colour = .dd_col_pop
      )
  }

  if (!is.null(curve_ind) && nrow(curve_ind)) {
    p <- p +
      ggplot2::geom_line(
        data = curve_ind,
        ggplot2::aes(x = .data$x, y = .data$.value, group = .data$.group),
        alpha = 0.3,
        linewidth = 0.4,
        colour = .dd_col_pop
      )
  }

  if (!is.null(curve_pop$.lower)) {
    if (multi_group) {
      p <- p +
        ggplot2::geom_ribbon(
          data = curve_pop,
          ggplot2::aes(
            x = .data$x,
            ymin = .data$.lower,
            ymax = .data$.upper,
            fill = .data$.group
          ),
          alpha = 0.18,
          colour = NA
        )
    } else {
      p <- p +
        ggplot2::geom_ribbon(
          data = curve_pop,
          ggplot2::aes(x = .data$x, ymin = .data$.lower, ymax = .data$.upper),
          alpha = 0.18,
          fill = .dd_col_accent,
          colour = NA
        )
    }
  }

  if (multi_group) {
    p <- p +
      ggplot2::geom_line(
        data = curve_pop,
        ggplot2::aes(x = .data$x, y = .data$.value, colour = .data$.group),
        linewidth = 1.1
      )
  } else {
    p <- p +
      ggplot2::geom_line(
        data = curve_pop,
        ggplot2::aes(x = .data$x, y = .data$.value, group = .data$.group),
        linewidth = 1.2,
        colour = .dd_col_pop
      )
  }

  p +
    .dd_x_scale(x_trans) +
    .dd_plot_theme() +
    ggplot2::labs(x = xlab, y = ylab, title = title, colour = NULL, fill = NULL)
}

# Apply a user `which` filter to the random-effect term names, erroring if the
# filter matches nothing (otherwise pivot_longer() fails opaquely downstream).
.dd_filter_qq_terms <- function(recols, which) {
  if (is.null(which)) {
    return(recols)
  }
  kept <- intersect(recols, which)
  if (!length(kept)) {
    cli::cli_abort(c(
      "{.arg which} matched no random-effect term.",
      "i" = "Available term{?s}: {.val {recols}}."
    ))
  }
  kept
}

# Random-effect QQ from a long (term, value) frame. Shared by the tmb/choice
# plot_qq methods.
.dd_qq_plot <- function(long, title = NULL, min_n = 3L) {
  n_per <- tapply(long$.value, long$.term, function(v) sum(is.finite(v)))
  if (any(n_per < min_n)) {
    cli::cli_warn(c(
      "Fewer than {min_n} subjects for a random-effect QQ plot.",
      "i" = "The plot is shown but the normal reference is unreliable with so few subjects."
    ))
  }
  ggplot2::ggplot(long, ggplot2::aes(sample = .data$.value)) +
    ggplot2::stat_qq(alpha = 0.5, colour = .dd_col_pop) +
    ggplot2::stat_qq_line(linetype = 2, colour = .dd_col_accent) +
    ggplot2::facet_wrap(ggplot2::vars(.data$.term), scales = "free") +
    .dd_plot_theme() +
    ggplot2::labs(
      x = "Theoretical quantiles",
      y = "Sample quantiles (random effects)",
      title = title
    )
}

# ==============================================================================
# TMB mixed-effects indifference-point model
# ==============================================================================

#' Plot a mixed-effects discounting model
#'
#' Visualize a fitted [fit_dd_tmb()] model. `type = "population"` draws the
#' population (random-effects-at-zero) discount curve over the observed
#' indifference points; `type = "individual"` adds the per-subject curves (the
#' shrinkage picture). For a fit with factors the population curve is drawn once
#' per factor-level combination; use `at` to condition on specific levels or
#' covariate values. `type = "parameters"` shows the distribution of the
#' subject-specific discount rate `k` (log scale); `type = "resid"` plots
#' standardized (Pearson) residuals against fitted values.
#'
#' @param x A `beezdiscounting_tmb` object.
#' @param type One of `"population"`, `"individual"`, `"parameters"`, `"resid"`.
#' @param ids Optional subset of subject ids for `type = "individual"`.
#' @param at Optional named list conditioning the population curve on factor
#'   levels / covariate values (e.g. `list(group = "A")`), passed to the same
#'   reference-grid machinery as [get_dd_comparisons()].
#' @param n_points Number of delay points in the curve grid.
#' @param x_trans Delay-axis scale: `"log10"` (default) or `"linear"`.
#' @param show_observed Overlay the observed indifference points.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' \donttest{
#' sim <- simulate_dd_ip(n_subjects = 12, seed = 1)
#' fit <- fit_dd_tmb(sim, equation = "mazur")
#' plot(fit)
#' plot(fit, type = "individual")
#' }
#' @export
plot.beezdiscounting_tmb <- function(
  x,
  type = c("population", "individual", "parameters", "resid"),
  ids = NULL,
  at = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
) {
  type <- match.arg(type)
  x_trans <- match.arg(x_trans)
  fit <- x

  if (type == "parameters") {
    return(.dd_plot_k_distribution(fit))
  }
  if (type == "resid") {
    return(.dd_plot_resid(fit))
  }

  .dd_ip_curve_plot(fit, type, ids, at, n_points, x_trans, show_observed)
}

# Shared indifference-point curve plot (population +/- per-subject + observed),
# used by both the TMB and brms IP tiers. The brms tier's population curve
# carries a credible band via .dd_ip_curve_values().
.dd_ip_curve_plot <- function(
  fit,
  type,
  ids,
  at,
  n_points,
  x_trans,
  show_observed
) {
  x_var <- .dd_x_var(fit)
  y_var <- fit$param_info$y_var
  nd_pop <- .dd_curve_newdata(
    fit,
    at = at,
    n_points = n_points,
    x_trans = x_trans
  )
  curve_pop <- .dd_ip_curve_values(fit, nd_pop, level = "population")

  curve_ind <- NULL
  if (type == "individual") {
    nd_ind <- .dd_individual_newdata(
      fit,
      ids = ids,
      n_points = n_points,
      x_trans = x_trans
    )
    curve_ind <- .dd_ip_curve_values(fit, nd_ind, level = "subject")
  }

  observed <- NULL
  if (show_observed) {
    observed <- data.frame(x = fit$data[[x_var]], y = fit$data[[y_var]])
    observed <- .dd_drop_nonpos_x(observed, "x", x_trans, "indifference points")
  }

  .dd_assemble_curve(
    curve_pop,
    curve_ind,
    observed,
    x_trans = x_trans,
    xlab = "Delay",
    ylab = "Indifference point"
  )
}

# Subject-k distribution (log10 axis: k is right-skewed on the natural scale).
.dd_plot_k_distribution <- function(fit) {
  sp <- stats::predict(fit, type = "parameters")
  ggplot2::ggplot(sp, ggplot2::aes(x = .data$k)) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = .dd_col_pop,
      colour = "white",
      alpha = 0.9
    ) +
    ggplot2::scale_x_log10() +
    .dd_plot_theme() +
    ggplot2::labs(x = "Subject discount rate k (log scale)", y = "Count")
}

# Standardized (Pearson) residuals vs fitted. Exploratory: for the bounded
# SLT-beta likelihood these are not expected to be normal, but systematic trend
# or funneling still flags lack of fit.
.dd_plot_resid <- function(fit) {
  aug <- generics::augment(fit)
  ggplot2::ggplot(aug, ggplot2::aes(x = .data$.fitted, y = .data$.std_resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = .dd_col_accent) +
    ggplot2::geom_point(alpha = 0.3, colour = .dd_col_pop) +
    .dd_plot_theme() +
    ggplot2::labs(
      x = "Fitted indifference point",
      y = "Standardized (Pearson) residual"
    )
}

#' Random-effect normal QQ plots for discounting models
#'
#' QQ plot of the estimated subject random-effect deviates against a normal
#' reference -- the standard check on the Gaussian random-effects assumption.
#' Methods exist for the TMB indifference-point ([fit_dd_tmb()]) and choice
#' ([fit_dd_choice()]) models. Bayesian fits are intentionally excluded; use
#' `brms::pp_check()` and MCMC diagnostics there instead. The plotted deviates
#' are shrunken (empirical-Bayes) estimates, so the normal reference is
#' approximate.
#'
#' @param object A fitted `beezdiscounting_tmb` or `beezdiscounting_choice`
#'   model.
#' @param which Optional character vector selecting random-effect terms to show.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @name plot_qq.beezdiscounting
NULL

#' @rdname plot_qq.beezdiscounting
#' @exportS3Method beezdemand::plot_qq
plot_qq.beezdiscounting_tmb <- function(object, which = NULL, ...) {
  re <- nlme::ranef(object)
  recols <- intersect(c("u_i", "re_k", "re_phi", "re_s"), names(re))
  if (!length(recols)) {
    cli::cli_abort("No random-effect deviates available for a QQ plot.")
  }
  recols <- .dd_filter_qq_terms(recols, which)
  long <- tidyr::pivot_longer(
    re[c("id", recols)],
    cols = dplyr::all_of(recols),
    names_to = ".term",
    values_to = ".value"
  )
  .dd_qq_plot(long, title = "Random-effect normal QQ")
}

# ==============================================================================
# Choice (SS-vs-LL) model -- structural and descriptive
# ==============================================================================

# Structural implied discount curve: per-row k (population or subject) mapped
# through the discount function. Choice fits are 1-parameter (mazur/exponential)
# so s = 1.
.dd_choice_curve_values <- function(fit, nd, level) {
  x_var <- .dd_x_var(fit)
  k <- .dd_tmb_predict_k(
    fit,
    nd[setdiff(names(nd), ".group")],
    level = if (level == "population") "population" else "subject"
  )
  data.frame(
    x = nd[[x_var]],
    .value = .dd_discount_mu(k, nd[[x_var]], fit$param_info$equation),
    .group = nd$.group
  )
}

# Calibration: observed choice proportion vs fitted P(LL), binned by fitted
# probability. Well-defined for structural and descriptive fits regardless of
# the per-trial amount/delay geometry.
.dd_plot_calibration <- function(fit, bins = 10L) {
  # TMB choice predict() returns .prob; the brms choice predict() returns the
  # posterior-median P(LL) in .fitted via type = "response".
  pr <- if (inherits(fit, "beezdiscounting_choice_brms")) {
    stats::predict(fit, type = "response", level = "subject")$.fitted
  } else {
    stats::predict(fit, type = "prob", level = "subject")$.prob
  }
  obs <- fit$data$choice
  brks <- seq(0, 1, length.out = bins + 1L)
  binf <- cut(pr, breaks = brks, include.lowest = TRUE)
  df <- data.frame(pred = pr, observed = obs, binf = binf)
  agg <- do.call(
    rbind,
    lapply(split(df, df$binf), function(d) {
      if (!nrow(d)) {
        return(NULL)
      }
      data.frame(pred = mean(d$pred), observed = mean(d$observed), n = nrow(d))
    })
  )
  ggplot2::ggplot(agg, ggplot2::aes(x = .data$pred, y = .data$observed)) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      colour = .dd_col_accent
    ) +
    ggplot2::geom_line(colour = .dd_col_pop, alpha = 0.6) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$n),
      colour = .dd_col_pop,
      alpha = 0.8
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    .dd_plot_theme() +
    ggplot2::labs(
      x = "Predicted P(choose LL)",
      y = "Observed proportion (binned)",
      size = "Trials"
    )
}

#' Plot a choice-based discounting model
#'
#' Visualize a fitted [fit_dd_choice()] model. For a *structural* fit,
#' `type = "population"` (default) draws the implied discount curve from the
#' estimated discount rate, `type = "individual"` adds the per-subject curves,
#' and `type = "parameters"` shows the subject-`k` distribution. For both
#' structural and *descriptive* (Young 2018) fits, `type = "calibration"` plots
#' the observed choice proportion against the fitted P(choose LL), binned by
#' fitted probability (a likelihood-geometry-agnostic goodness-of-fit view).
#' Descriptive fits have no structural discount rate, so the curve/parameter
#' types are unavailable.
#'
#' @param x A `beezdiscounting_choice` object.
#' @param type One of `"population"`, `"individual"`, `"calibration"`,
#'   `"parameters"`. Defaults to `"population"` (structural) or `"calibration"`
#'   (descriptive).
#' @param ids Optional subset of subject ids for `type = "individual"`.
#' @param at Optional named list conditioning the curve on factor levels /
#'   covariate values.
#' @param n_points Number of delay points in the curve grid.
#' @param x_trans Delay-axis scale: `"log10"` (default) or `"linear"`.
#' @param show_observed Unused for the implied discount curve (kept for a
#'   consistent signature across tiers).
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' \donttest{
#' ch <- simulate_dd_choice(n_subjects = 12, seed = 1)
#' fit <- fit_dd_choice(ch, equation = "mazur")
#' plot(fit)
#' plot(fit, type = "calibration")
#' }
#' @export
plot.beezdiscounting_choice <- function(
  x,
  type = NULL,
  ids = NULL,
  at = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
) {
  fit <- x
  x_trans <- match.arg(x_trans)
  mode <- fit$param_info$mode
  if (is.null(type)) {
    type <- if (identical(mode, "descriptive")) "calibration" else "population"
  }
  type <- match.arg(
    type,
    c("population", "individual", "calibration", "parameters")
  )

  if (
    identical(mode, "descriptive") &&
      type %in% c("population", "individual", "parameters")
  ) {
    cli::cli_abort(c(
      "{.val {type}} is not available for a descriptive choice fit.",
      "i" = "A descriptive (Young 2018) fit has no structural discount rate.",
      ">" = "Use {.code type = \"calibration\"} or {.fn plot_qq}."
    ))
  }

  if (type == "calibration") {
    return(.dd_plot_calibration(fit))
  }
  if (type == "parameters") {
    return(.dd_plot_k_distribution(fit))
  }

  nd_pop <- .dd_curve_newdata(
    fit,
    at = at,
    n_points = n_points,
    x_trans = x_trans
  )
  curve_pop <- .dd_choice_curve_values(fit, nd_pop, level = "population")
  curve_ind <- NULL
  if (type == "individual") {
    nd_ind <- .dd_individual_newdata(
      fit,
      ids = ids,
      n_points = n_points,
      x_trans = x_trans
    )
    curve_ind <- .dd_choice_curve_values(fit, nd_ind, level = "subject")
  }

  .dd_assemble_curve(
    curve_pop,
    curve_ind,
    observed = NULL,
    x_trans = x_trans,
    xlab = "Delay",
    ylab = "Relative subjective value"
  )
}

#' @rdname plot_qq.beezdiscounting
#' @exportS3Method beezdemand::plot_qq
plot_qq.beezdiscounting_choice <- function(object, which = NULL, ...) {
  mode <- object$param_info$mode
  re <- nlme::ranef(object)
  recols <- if (identical(mode, "descriptive")) {
    intersect(c("b_mag", "b_delay"), names(re))
  } else {
    intersect("u_i", names(re))
  }
  if (!length(recols)) {
    cli::cli_abort(c(
      "No subject random effects to plot.",
      "i" = "A pooled descriptive fit ({.code random_slopes = FALSE}) has no random effects."
    ))
  }
  recols <- .dd_filter_qq_terms(recols, which)
  long <- tidyr::pivot_longer(
    re[c("id", recols)],
    cols = dplyr::all_of(recols),
    names_to = ".term",
    values_to = ".value"
  )
  .dd_qq_plot(long, title = "Random-effect normal QQ")
}

# ==============================================================================
# Group comparison forest plot
# ==============================================================================

#' Plot group differences in discount rate
#'
#' Forest plot of the pairwise (or treatment-vs-control) contrasts from
#' [get_dd_comparisons()]. `type = "ratio"` (default) shows each k ratio on a
#' log axis with a reference line at 1; `type = "difference"` shows the
#' difference in log10 k with a reference line at 0. A contrast is flagged when
#' its interval excludes the null -- a backend-agnostic encoding that works for
#' both the frequentist (TMB) and Bayesian (brms) backends.
#'
#' @param x A `beezdiscounting_comparison` object.
#' @param type `"ratio"` (default) or `"difference"`.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' \donttest{
#' sim <- simulate_dd_ip(n_subjects = 30, n_conditions = 2, delta_k = c(0, 0.5), seed = 1)
#' fit <- fit_dd_tmb(sim, equation = "mazur", factors = "condition")
#' cmp <- get_dd_comparisons(fit)
#' plot(cmp)
#' }
#' @export
plot.beezdiscounting_comparison <- function(
  x,
  type = c("ratio", "difference"),
  ...
) {
  type <- match.arg(type)
  td <- generics::tidy(x, exponentiate = (type == "ratio"))
  if (!nrow(td)) {
    cli::cli_abort(c(
      "No contrasts to plot.",
      "i" = "The comparison object has no rows (was the fit intercept-only?)."
    ))
  }

  stat_cols <- c(
    "param",
    "contrast",
    "estimate",
    "std.error",
    "statistic",
    "df",
    "conf.low",
    "conf.high",
    "p.value",
    "post.prob"
  )
  by_cols <- setdiff(names(td), stat_cols)
  td$.label <- if (length(by_cols)) {
    paste0(
      do.call(
        paste,
        c(
          lapply(by_cols, function(cc) as.character(td[[cc]])),
          list(sep = ", ")
        )
      ),
      " | ",
      td$contrast
    )
  } else {
    td$contrast
  }
  td$.label <- factor(td$.label, levels = rev(unique(td$.label)))

  null_val <- if (type == "ratio") 1 else 0
  excl <- (td$conf.low > null_val & td$conf.high > null_val) |
    (td$conf.low < null_val & td$conf.high < null_val)
  td$.flag <- factor(
    ifelse(excl, "CI excludes null", "CI includes null"),
    levels = c("CI excludes null", "CI includes null")
  )

  p <- ggplot2::ggplot(
    td,
    ggplot2::aes(x = .data$estimate, y = .data$.label, colour = .data$.flag)
  ) +
    ggplot2::geom_vline(
      xintercept = null_val,
      linetype = 2,
      colour = "grey50"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high)
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(
      values = c(
        "CI excludes null" = .dd_col_accent,
        "CI includes null" = "grey55"
      ),
      drop = FALSE
    ) +
    .dd_plot_theme() +
    ggplot2::labs(
      x = if (type == "ratio") {
        "k ratio (log scale)"
      } else {
        "Difference in log10 k"
      },
      y = NULL,
      colour = NULL
    )
  if (type == "ratio") {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}

# ==============================================================================
# Bayesian (brms) tiers
# ==============================================================================

# Posterior subject-k caterpillar (point + credible interval per subject) on a
# log10 axis. Shared by both brms tiers.
.dd_plot_k_caterpillar <- function(fit) {
  sp <- stats::predict(fit, type = "parameters")
  if (is.null(sp)) {
    cli::cli_abort(c(
      "The subject-{.field k} plot is unavailable for this fit.",
      "i" = "Subject-level k is undefined when the fixed-effect design varies within subject."
    ))
  }
  sp <- sp[order(sp$k), , drop = FALSE]
  sp$id <- factor(as.character(sp$id), levels = as.character(sp$id))
  p <- ggplot2::ggplot(sp, ggplot2::aes(y = .data$id, x = .data$k)) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$k_lower, xmax = .data$k_upper),
      colour = .dd_col_pop,
      alpha = 0.6
    ) +
    ggplot2::geom_point(colour = .dd_col_pop) +
    ggplot2::scale_x_log10() +
    .dd_plot_theme() +
    ggplot2::labs(x = "Subject discount rate k (log scale)", y = "Subject")
  # Per-subject id labels become an illegible smear once there are many
  # subjects; drop them and label the axis as a rank instead.
  if (nrow(sp) > 30L) {
    p <- p +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(y = "Subject (ordered by k)")
  }
  p
}

#' Plot a Bayesian (brms) indifference-point discounting model
#'
#' Visualize a fitted [fit_dd_brms()] model. `type = "population"` (default)
#' draws the posterior-median discount curve with a credible band for the
#' population conditional mean (random effects set to zero -- a credible band,
#' not a posterior-predictive band) over the observed (boundary-squeezed)
#' indifference points; `type = "individual"` adds per-subject median curves;
#' `type = "parameters"` shows the posterior subject-`k` caterpillar;
#' `type = "resid"` plots standardized residuals against fitted values.
#'
#' @param x A `beezdiscounting_brms` object.
#' @param type One of `"population"`, `"individual"`, `"parameters"`, `"resid"`.
#' @param ids Optional subset of subject ids for `type = "individual"`.
#' @param at Optional named list conditioning the population curve on factor
#'   levels / covariate values.
#' @param n_points Number of delay points in the curve grid.
#' @param x_trans Delay-axis scale: `"log10"` (default) or `"linear"`.
#' @param show_observed Overlay the observed indifference points.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @export
plot.beezdiscounting_brms <- function(
  x,
  type = c("population", "individual", "parameters", "resid"),
  ids = NULL,
  at = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
) {
  fit <- x
  type <- match.arg(type)
  x_trans <- match.arg(x_trans)
  if (type == "parameters") {
    return(.dd_plot_k_caterpillar(fit))
  }
  if (type == "resid") {
    return(.dd_plot_resid(fit))
  }
  .dd_ip_curve_plot(fit, type, ids, at, n_points, x_trans, show_observed)
}

# Design matrix for a brms fit's k-model on newdata (factors releveled to the
# fitted levels), aligned to the fitted columns. Used to map posterior log-k
# draws to an implied discount curve.
.dd_brms_design_matrix <- function(fit, newdata) {
  rhs <- fit$formula_details$rhs %||% stats::as.formula("~ 1")
  ca <- fit$formula_details$contrasts
  for (f in fit$param_info$factors) {
    if (
      !is.null(f) &&
        nzchar(f) &&
        f %in% names(newdata) &&
        f %in% names(fit$data)
    ) {
      newdata[[f]] <- factor(
        newdata[[f]],
        levels = levels(as.factor(fit$data[[f]]))
      )
    }
  }
  X <- stats::model.matrix(rhs, data = newdata, contrasts.arg = ca)
  X[, colnames(fit$formula_details$X), drop = FALSE]
}

# Posterior implied discount curve (population, random effects at zero) for a
# Bayesian choice fit: map population log-k draws through the discount function
# over the delay grid and summarize to median + credible band.
.dd_choice_brms_population_curve <- function(fit, nd, probs = c(0.025, 0.975)) {
  draws <- .dd_brms_draws_matrix(fit)
  kvars <- .dd_brms_logk_draw_vars(fit, colnames(draws))
  b <- as.matrix(draws[, kvars, drop = FALSE])
  X <- .dd_brms_design_matrix(fit, nd[setdiff(names(nd), ".group")])
  k <- exp(b %*% t(X)) # [ndraws x nrow(nd)]
  xv <- nd[[.dd_x_var(fit)]]
  xmat <- matrix(xv, nrow = nrow(k), ncol = length(xv), byrow = TRUE)
  V <- .dd_discount_mu(k, xmat, fit$param_info$equation)
  data.frame(
    x = xv,
    .value = apply(V, 2, stats::median),
    .lower = apply(V, 2, stats::quantile, probs = probs[1]),
    .upper = apply(V, 2, stats::quantile, probs = probs[2]),
    .group = nd$.group
  )
}

# Per-subject implied discount curves (posterior-median k) for a Bayesian
# choice fit.
.dd_choice_brms_individual <- function(fit, ids, n_points, x_trans) {
  x_var <- .dd_x_var(fit)
  xs <- .dd_x_seq(fit$data[[x_var]], n_points, x_trans)
  sp <- fit$subject_pars
  if (is.null(sp)) {
    cli::cli_abort(c(
      "Per-subject curves are unavailable for this fit.",
      "i" = "Subject-level k is undefined when the fixed-effect design varies within subject."
    ))
  }
  all_ids <- as.character(sp$id)
  if (!is.null(ids)) {
    ids <- unique(as.character(ids))
    if (!length(ids)) {
      cli::cli_abort("{.arg ids} is empty.")
    }
    unknown <- setdiff(ids, all_ids)
    if (length(unknown)) {
      cli::cli_abort("Unknown {.arg ids}: {.val {unknown}}.")
    }
    sp <- sp[match(ids, all_ids), , drop = FALSE]
    all_ids <- ids
  }
  rows <- lapply(seq_len(nrow(sp)), function(i) {
    data.frame(
      x = xs,
      .value = .dd_discount_mu(sp$k[i], xs, fit$param_info$equation),
      .group = sp$id[i]
    )
  })
  out <- do.call(rbind, rows)
  out$.group <- factor(as.character(out$.group), levels = all_ids)
  out
}

#' Plot a Bayesian (brms) choice-based discounting model
#'
#' Visualize a fitted [fit_dd_choice_brms()] model. `type = "population"`
#' (default) draws the posterior implied discount curve with a credible band;
#' `type = "individual"` adds per-subject (posterior-median) curves;
#' `type = "parameters"` shows the posterior subject-`k` caterpillar; and
#' `type = "calibration"` plots observed choice proportion against fitted
#' P(choose LL), binned by fitted probability.
#'
#' @param x A `beezdiscounting_choice_brms` object.
#' @param type One of `"population"`, `"individual"`, `"calibration"`,
#'   `"parameters"`.
#' @param ids Optional subset of subject ids for `type = "individual"`.
#' @param at Optional named list conditioning the population curve on factor
#'   levels / covariate values.
#' @param n_points Number of delay points in the curve grid.
#' @param x_trans Delay-axis scale: `"log10"` (default) or `"linear"`.
#' @param show_observed Unused for the implied discount curve.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @export
plot.beezdiscounting_choice_brms <- function(
  x,
  type = c("population", "individual", "calibration", "parameters"),
  ids = NULL,
  at = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
) {
  fit <- x
  type <- match.arg(type)
  x_trans <- match.arg(x_trans)
  if (type == "calibration") {
    return(.dd_plot_calibration(fit))
  }
  if (type == "parameters") {
    return(.dd_plot_k_caterpillar(fit))
  }
  nd_pop <- .dd_curve_newdata(
    fit,
    at = at,
    n_points = n_points,
    x_trans = x_trans
  )
  curve_pop <- .dd_choice_brms_population_curve(fit, nd_pop)
  curve_ind <- if (type == "individual") {
    .dd_choice_brms_individual(fit, ids, n_points, x_trans)
  } else {
    NULL
  }
  .dd_assemble_curve(
    curve_pop,
    curve_ind,
    observed = NULL,
    x_trans = x_trans,
    xlab = "Delay",
    ylab = "Relative subjective value"
  )
}
