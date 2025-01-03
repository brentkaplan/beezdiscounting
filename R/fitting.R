
#' Fit Delay-Discounting Model
#'
#' This function fits a delay-discounting model to the given dataset using the specified equation and method.
#'
#' @param dat A data frame containing delay (`x`) and indifference point (`y`) data. For "two stage" methods, the data must include an `id` column to identify participants.
#' @param equation A character string specifying the delay-discounting equation to use. Options include:
#'   - `"mazur"` or `"hyperbolic"`: Hyperbolic delay-discounting model (\eqn{y = 1 / (1 + k \cdot x)}).
#'   - `"exponential"`: Exponential delay-discounting model (\eqn{y = \exp(-k \cdot x)}).
#' @param method A character string specifying the method for fitting the model. Options include:
#'   - `"pooled"` or `"agg"`: Fits the model using pooled data.
#'   - `"mean"`: Fits the model using the mean of indifference points at each delay.
#'   - `"ts"` or `"two stage"`: Fits the model separately for each participant (requires an `id` column in `dat`).
#'
#' @return A list object of class `"fit_dd"`, containing:
#'   - The fitted model(s).
#'   - The original dataset (`dat`).
#'   - The specified method (`method`).
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:2, each = 6),
#'   x = rep(c(1, 7, 30, 90, 180, 365), 2),
#'   y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
#' )
#' fit_dd(data, equation = "mazur", method = "two stage")
fit_dd <- function(
  dat,
  equation,
  method
) {

  stopifnot(any(
    equation %in% c(
      "mazur",
      "Mazur",
      "hyperbolic",
      "Hyperbolic",
      "exponential",
      "Exponential"
    )
  ))

  stopifnot(any(
    method %in% c(
      "pooled",
      "Pooled",
      "agg",
      "Agg",
      "mean",
      "Mean",
      "ts",
      "two stage",
      "Two Stage"
    )
    ))

  dat$id <- factor(dat$id)

  if (any(equation %in% c("mazur", "Mazur", "hyperbolic", "Hyperbolic"))) {
    fo <- y ~ 1 / (1 + k * x)
  } else if (any(equation %in% c("exponential", "Exponential"))) {
    fo <- y ~ exp(-k * x)
  }

  if (any(method %in% c("pooled", "Pooled", "agg", "Agg"))) {
    fit <- minpack.lm::nlsLM(
      formula = fo,
      start = list(k = 0.01),
      data = dat
    )
  } else if (any(method %in% c("mean", "Mean"))) {

    dat <- aggregate(y ~ x, data = dat, mean)
    fit <- minpack.lm::nlsLM(
      formula = fo,
      start = list(k = 0.01),
      data = dat
    )

  } else if (any(method %in% c("ts", "two stage", "Two Stage"))) {

   fit_model <- function(data, formula, start_params) {
    minpack.lm::nlsLM(
      formula = formula,
      start = start_params,
      data = data
    )
  }

    safe_fit_model <- purrr::safely(fit_model)

    data_split <- split(dat, dat$id)

    fit <- data_split |>
      purrr::map(~ safe_fit_model(.x, formula = fo, start_params = list(k = 0.01)))
  }

  dd_out <- list(
      fit,
      dat,
      method
    )
  class(dd_out) <- c("fit_dd", class(dd_out))

  return(dd_out)

}


#' Plot Delay-Discounting Model
#'
#' This function generates a plot of the delay-discounting data and the fitted model.
#'
#' @param fit_dd_object A fitted delay-discounting model object of class `"fit_dd"`, created by the `fit_dd()` function.
#' @param xlabel A character string specifying the label for the x-axis. Default is `"Delay"`.
#' @param ylabel A character string specifying the label for the y-axis. Default is `"Indifference Point"`.
#' @param title A character string specifying the plot title. Default is `""`.
#' @param logx Logical. If `TRUE`, the x-axis is log-transformed. Default is `TRUE`.
#'
#' @return A ggplot object representing the fitted model and data.
#' @importFrom stats aggregate predict
#' @export
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:2, each = 6),
#'   x = rep(c(1, 7, 30, 90, 180, 365), 2),
#'   y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
#' )
#' fit <- fit_dd(data, equation = "mazur", method = "mean")
#' plot_dd(fit)
plot_dd <- function(
  fit_dd_object,
  xlabel = "Delay",
  ylabel = "Indifference Point",
  title = "",
  logx = TRUE
) {

  stopifnot(any(class(fit_dd_object) %in% "fit_dd"))

  new_x <- seq(min(fit_dd_object[[2]]$x), max(fit_dd_object[[2]]$x), length.out = 100)

  if (fit_dd_object[[3]] %in% c("pooled", "Pooled", "agg", "Agg")) {
    pred <- predict(fit_dd_object[[1]], newdata = data.frame(x = new_x))
    plt <- ggplot2::ggplot(
      fit_dd_object[[2]],
      ggplot2::aes(x = x, y = y)
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_line(
        data = data.frame(x = new_x, y = pred),
        ggplot2::aes(x = x, y = y),
        color = "red"
      )


  } else if (fit_dd_object[[3]] %in% c("mean", "Mean")) {
    pred <- predict(fit_dd_object[[1]], newdata = data.frame(x = new_x))
    dat <- aggregate(y ~ x, data = fit_dd_object[[2]], mean)
    plt <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = x, y = y)
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_line(
        data = data.frame(x = new_x, y = pred),
        ggplot2::aes(x = x, y = y),
        color = "red"
      )

  } else if (fit_dd_object[[3]] %in% c("ts", "two stage", "Two Stage")) {

    pred <- purrr::map2_dfr(
      .x = fit_dd_object[[1]],
      .y = names(fit_dd_object[[1]]),
      ~ {
        model <- .x$result
        if (is.null(model)) {
          return(NULL)  # Skip if the model is NULL
        }
        preds <- predict(model, newdata = data.frame(x = new_x))
        tibble::tibble(
          id = .y,
          new_x = new_x,
          pred = preds
        )
      }
    )

    plt <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = pred,
        alpha = .5,
        ggplot2::aes(x = new_x, y = pred, group = id)
      ) +
      ggplot2::theme(
        legend.position = "none"
      )
  }

  plt <- plt +
    beezdemand::theme_apa() +
    ggplot2::labs(
      x = xlabel,
      y = ylabel,
      title = title
    )

  if (logx) {
    plt <- plt +
      ggplot2::scale_x_log10()
  }

  return(plt)
}


#' Extract Results from Delay-Discounting Model
#'
#' This function extracts model parameter estimates, fit statistics, and confidence intervals from a fitted delay-discounting model.
#'
#' @param fit_dd_object A fitted delay-discounting model object of class `"fit_dd"`, created by the `fit_dd()` function.
#'
#' @return A tibble containing the following columns:
#'   - `id`: The participant or group ID (if applicable).
#'   - `term`: The model parameter (e.g., `k`).
#'   - `estimate`: The estimated value of the parameter.
#'   - `std.error`: The standard error of the parameter estimate.
#'   - `statistic`: The t-statistic for the parameter estimate.
#'   - `p.value`: The p-value for the parameter estimate.
#'   - `conf_low`: The lower bound of the 95% confidence interval.
#'   - `conf_high`: The upper bound of the 95% confidence interval.
#'   - `R2`: The coefficient of determination (\eqn{R^2}).
#' @export
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:2, each = 6),
#'   x = rep(c(1, 7, 30, 90, 180, 365), 2),
#'   y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.85, 0.55, 0.35, 0.15, 0.1, 0.05)
#' )
#' fit <- fit_dd(data, equation = "mazur", method = "two stage")
#' results_dd(fit)
results_dd <- function(fit_dd_object) {

  stopifnot(any(class(fit_dd_object) %in% "fit_dd"))

  if (fit_dd_object[[3]] %in% c("pooled", "Pooled", "agg", "Agg", "mean", "Mean")) {
    out <- broom::tidy(fit_dd_object[[1]]) |>
      dplyr::bind_cols(broom::glance(fit_dd_object[[1]])) |>
      dplyr::mutate(
        method = fit_dd_object[[3]],
        R2 = calc_r2(fit_dd_object[[1]]),
        model = fit_dd_object[1],
        conf = purrr::pmap(
          list(estimate = estimate, std_error = std.error, model = model),
          ~ calc_conf_int(..1, ..2, ..3, alpha = 0.05)
        ),
        conf_low = purrr::map_dbl(conf, 1),
        conf_high = purrr::map_dbl(conf, 2)
      ) |>
      dplyr::relocate(method, .before = term) |>
      dplyr::select(-conf, -model)

    if (fit_dd_object[[3]] %in% c("mean", "Mean")) {
      out <- out |>
        dplyr::mutate(calc_aucs(fit_dd_object[[2]])) |>
        dplyr::relocate(dplyr::starts_with("auc"), .after = R2)
    }

  } else if (fit_dd_object[[3]] %in% c("ts", "two stage", "Two Stage")) {
    fit_results <- tibble::tibble(
      id = names(fit_dd_object[[1]]),
      fit = purrr::map(fit_dd_object[[1]], "result"),
      error = purrr::map(fit_dd_object[[1]], "error")
    ) |>
      dplyr::mutate(
        tidy_summary = purrr::map(fit, ~ if (!is.null(.)) broom::tidy(.) else NA),
        R2 = purrr::map_dbl(fit, calc_r2),
        conf_int = purrr::map(fit, ~ {
          if (!is.null(.)) {
            tidy_res <- broom::tidy(.)
            purrr::pmap(
              list(
                estimate = tidy_res$estimate,
                std_error = tidy_res$std.error,
                model = list(.)
              ),
              ~ calc_conf_int(..1, ..2, ..3)
            )
          } else {
            NA
          }
        }),
        glance_summary = purrr::map(fit, ~ if (!is.null(.)) broom::glance(.) else NA)
      )

    fit_results <- fit_dd_object[[2]] |>
      dplyr::mutate(id = as.character(id)) |>
      dplyr::group_by(id) |>
      dplyr::group_split() |>
      purrr::map_dfr(calc_aucs) |>
      dplyr::left_join(x = fit_results, y = _, by = "id")

    out <- fit_results |>
      dplyr::select(id, tidy_summary, glance_summary, R2, tidyr::starts_with("auc"), conf_int) |>
      tidyr::unnest(cols = tidy_summary) |>
      tidyr::unnest(cols = glance_summary) |>
      dplyr::mutate(
        conf_low = purrr::map_dbl(conf_int, ~ .[[1]][1]),
        conf_high = purrr::map_dbl(conf_int, ~ .[[1]][2])
      ) |>
      dplyr::select(-conf_int) |>
      dplyr::mutate(method = fit_dd_object[[3]]) |>
      dplyr::relocate(method, .before = id)


  }

  return(out)

}
