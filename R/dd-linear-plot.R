# ==============================================================================
# plot() for the linearized Mazur tier (beezdiscounting_linear). Lives in its
# own file rather than dd-plot.R so the closed-form tier's plotting surface
# can be read in one place; it reuses dd-plot.R's shared assembler, palette,
# theme, caterpillar, and residual helpers so the pictures match the TMB and
# brms tiers.
# ==============================================================================

#' Plot a linearized Mazur discounting fit
#'
#' Visualize a fitted [fit_dd_linear()] model. `type = "population"` (default)
#' draws the hyperbola implied by each condition's geometric-mean discount rate
#' `exp(mu)` over the observed indifference points (one curve per condition when
#' the fit has a factor); `type = "individual"` adds each subject's own
#' closed-form hyperbola `1 / (1 + k_i x)` as thin lines; `type =
#' "transformed"` shows the linearization itself, `ln(1/D - 1)` against
#' `ln(delay)` per subject with a slope-1 line at that subject's `ln k`;
#' `type = "parameters"` shows the subject discount rates `k` with their
#' t-intervals (by condition, with the condition's geometric-mean `k` and its
#' interval overlaid, when the fit has a factor; otherwise ordered as a
#' caterpillar); `type = "resid"` plots the standardized residual on the
#' transformed scale against the raw-scale fitted indifference point.
#'
#' @details
#' The arguments match [plot.beezdiscounting_tmb()] except that `at` is absent:
#' a linearized fit has no covariates or reference grid, so there is nothing to
#' condition on. `x_trans`, `n_points` and `show_observed` apply to the
#' `"population"` and `"individual"` curves only; `ids` applies to
#' `"individual"` (default: every subject) and `"transformed"` (default: the
#' first 12 subjects, with a message when the fit has more).
#'
#' Unlike the mixed-model tiers there is no shrinkage: every subject's curve is
#' its own closed-form estimate, and the per-subject pictures do not need the
#' random-effects component. `"population"` requires it (the condition means
#' `mu`); when the design was unbalanced and `re` is `NULL`, `"population"`
#' errors, `"individual"` draws the subject curves only, and `"parameters"`
#' omits the condition means with a message.
#'
#' The `"transformed"` abscissa is `ln(delay)` by construction (a `log10` axis
#' would break the slope-1 reference), so `x_trans` is ignored there.
#' Indifference points at exactly 0 or 1 are shown at their observed value on
#' the raw scale, and at their `boundary`-handled value `d_used` on the
#' transformed scale; points dropped under `boundary = "drop"` do not appear on
#' the transformed scale. In `"parameters"`, `k` is drawn on a log10 axis (the
#' package's convention for subject discount rates); `coef()` and `confint()`
#' report the same quantities on the natural-log scale. `"resid"` uses
#' `augment()`'s `.std_resid`, the transformed-scale residual divided by the
#' model's error standard deviation; the abscissa is the raw-scale `.fitted`,
#' which is monotone in delay within a subject, so delay-dependent misfit shows
#' as a trend.
#'
#' @param x A `beezdiscounting_linear` object.
#' @param type One of `"population"`, `"individual"`, `"transformed"`,
#'   `"parameters"`, `"resid"`.
#' @param ids Optional subset of subject ids for `type = "individual"` and
#'   `type = "transformed"`.
#' @param n_points Number of delay points in the curve grid.
#' @param x_trans Delay-axis scale: `"log10"` (default) or `"linear"`.
#' @param show_observed Overlay the observed indifference points.
#' @param ... Unused.
#' @return A [ggplot2::ggplot] object.
#' @seealso [fit_dd_linear()], [beezdiscounting_linear-methods],
#'   `vignette("linearized-mazur")`.
#' @examples
#' sim <- simulate_dd_linear(
#'   n_subjects = 10, delays = c(7, 30, 90, 365),
#'   mu = c(A = -6, B = -4.5), sigma2 = 2, g = 8, seed = 1
#' )
#' fit <- fit_dd_linear(sim, factors = "condition")
#' plot(fit)
#' plot(fit, type = "individual")
#' plot(fit, type = "transformed", ids = c("A_1", "B_1"))
#' plot(fit, type = "parameters")
#' @export
plot.beezdiscounting_linear <- function(
  x,
  type = c("population", "individual", "transformed", "parameters", "resid"),
  ids = NULL,
  n_points = 200,
  x_trans = c("log10", "linear"),
  show_observed = TRUE,
  ...
) {
  type <- match.arg(type)
  x_trans <- match.arg(x_trans)
  fit <- x

  switch(
    type,
    population = .dd_lin_plot_curve(
      fit,
      type,
      ids,
      n_points,
      x_trans,
      show_observed
    ),
    individual = .dd_lin_plot_curve(
      fit,
      type,
      ids,
      n_points,
      x_trans,
      show_observed
    ),
    transformed = .dd_lin_plot_transformed(fit, ids),
    parameters = .dd_lin_plot_parameters(fit),
    resid = .dd_plot_resid(fit)
  )
}

# Raw-scale curves through the shared assembler: population hyperbolae from
# exp(mu) per condition (+/- per-subject hyperbolae) over the observed points.
.dd_lin_plot_curve <- function(
  fit,
  type,
  ids,
  n_points,
  x_trans,
  show_observed
) {
  if (type == "population" && is.null(fit$re)) {
    cli::cli_abort(c(
      "Random-effects component was not fitted.",
      "i" = "Use {.code type = \"individual\"} or {.code type = \"transformed\"} for the per-subject fits."
    ))
  }
  xs <- .dd_x_seq(fit$data$x, n_points, x_trans)
  curve_pop <- .dd_lin_curve_pop(fit, xs)

  curve_ind <- NULL
  if (type == "individual") {
    curve_ind <- .dd_lin_curve_ind(fit, .dd_lin_resolve_ids(fit, ids), xs)
  }

  observed <- NULL
  if (show_observed) {
    observed <- data.frame(x = fit$data$x, y = fit$data$y)
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

# Population curve frame (x, .value, .group). One group per condition level
# when the fit has a factor, otherwise the single "Population" group (the
# label .dd_curve_newdata() uses). Zero rows when `re` is NULL so the
# assembler's population layer draws nothing.
.dd_lin_curve_pop <- function(fit, xs) {
  if (is.null(fit$re)) {
    return(data.frame(
      x = numeric(0),
      .value = numeric(0),
      .group = factor(character(0), levels = "Population")
    ))
  }
  mu <- fit$re$mu
  labels <- if (is.null(fit$design$factor)) "Population" else names(mu)
  data.frame(
    x = rep(xs, times = length(mu)),
    .value = as.vector(vapply(
      mu,
      function(m) 1 / (1 + exp(m) * xs),
      numeric(length(xs))
    )),
    .group = factor(rep(labels, each = length(xs)), levels = labels)
  )
}

# Per-subject curve frame for the given ids (already validated).
.dd_lin_curve_ind <- function(fit, ids, xs) {
  k <- fit$subjects$k[match(ids, as.character(fit$subjects$id))]
  data.frame(
    x = rep(xs, times = length(ids)),
    .value = as.vector(vapply(
      k,
      function(ki) 1 / (1 + ki * xs),
      numeric(length(xs))
    )),
    .group = factor(rep(ids, each = length(xs)), levels = ids)
  )
}

# `ids` validation with the same contract as .dd_individual_newdata(): NULL
# means every fitted subject; otherwise unique character ids, all of which
# must have been fitted.
.dd_lin_resolve_ids <- function(fit, ids) {
  all_ids <- as.character(fit$subjects$id)
  if (is.null(ids)) {
    return(all_ids)
  }
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
  ids
}

# The linearization picture: ln(1/D - 1) vs ln(delay) per subject with a
# slope-1 line at the subject's ln k. Faceted, so the default subject set is
# capped (with a message) rather than drawing hundreds of panels.
.dd_lin_plot_transformed <- function(fit, ids, max_default = 12L) {
  if (is.null(ids)) {
    ids <- as.character(fit$subjects$id)
    if (length(ids) > max_default) {
      cli::cli_inform(c(
        "i" = "Showing the first {max_default} of {length(ids)} subjects.",
        " " = "Use {.arg ids} to choose which subjects to draw."
      ))
      ids <- ids[seq_len(max_default)]
    }
  } else {
    ids <- .dd_lin_resolve_ids(fit, ids)
  }

  d <- fit$data
  d <- d[as.character(d$id) %in% ids & is.finite(d$y_lin), , drop = FALSE]
  pts <- data.frame(
    id = factor(as.character(d$id), levels = ids),
    x = log(d$x),
    y = log(1 / d$d_used - 1)
  )
  lines <- data.frame(
    id = factor(ids, levels = ids),
    logk = fit$subjects$logk[match(ids, as.character(fit$subjects$id))],
    slope = 1
  )

  ggplot2::ggplot() +
    ggplot2::geom_abline(
      data = lines,
      ggplot2::aes(intercept = .data$logk, slope = .data$slope),
      colour = .dd_col_pop,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data$x, y = .data$y),
      colour = .dd_col_pop,
      alpha = 0.6
    ) +
    ggplot2::facet_wrap(~id) +
    .dd_plot_theme() +
    ggplot2::labs(x = "ln(delay)", y = "ln(1/D - 1)")
}

# Subject k with t-intervals by condition, the condition's geometric-mean k
# (exp(mu)) and its interval overlaid; the caterpillar when there is no factor.
.dd_lin_plot_parameters <- function(fit) {
  if (is.null(fit$design$factor)) {
    return(.dd_plot_k_caterpillar(fit))
  }
  sp <- stats::predict(fit, type = "parameters")
  p <- ggplot2::ggplot(sp, ggplot2::aes(x = .data$condition, y = .data$k)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$k_lower, ymax = .data$k_upper),
      colour = .dd_col_pop,
      alpha = 0.5,
      size = 0.3,
      position = ggplot2::position_jitter(width = 0.2, height = 0, seed = 1)
    )

  if (is.null(fit$re)) {
    cli::cli_inform(c(
      "i" = "Random-effects component was not fitted; condition means are not drawn."
    ))
  } else {
    ci <- exp(stats::confint(fit, level = fit$conf_level))
    pop <- data.frame(
      condition = factor(names(fit$re$mu), levels = levels(sp$condition)),
      k = exp(unname(fit$re$mu)),
      k_lower = unname(ci[, 1]),
      k_upper = unname(ci[, 2])
    )
    p <- p +
      ggplot2::geom_pointrange(
        data = pop,
        ggplot2::aes(ymin = .data$k_lower, ymax = .data$k_upper),
        colour = .dd_col_accent,
        size = 0.8
      )
  }

  p +
    ggplot2::scale_y_log10() +
    .dd_plot_theme() +
    ggplot2::labs(x = "Condition", y = "Subject discount rate k (log scale)")
}
