# R/dd-validate.R

#' Validate and coerce IP-family long data to proportions in [0, 1]
#'
#' The single coercion choke-point for indifference-point (IP) mixed-effects
#' discounting. It remaps caller column names to the canonical `id, x, y` long
#' format (matching `fit_dd()`) **while retaining any factor/covariate columns
#' named in `extra_cols`** so a single model frame can be complete-cased once
#' downstream (row-coherence; see `.dd_tmb_prepare_data`). It detects and
#' divides percent- or amount-scaled responses to the `[0, 1]` proportion
#' scale, and clamps mild post-coercion overshoot. **All coercion is loud**:
#' percent/amount division and any clamping emit a `warning()` that names how
#' many values were affected. This guarantees the downstream likelihood never
#' sees out-of-range `y` (we filter and coerce here rather than masking inside
#' the likelihood; see the design spec, "Implementation landmines").
#'
#' Scale detection:
#' \itemize{
#'   \item `response_scale = "proportion"` (default): if `max(y, na.rm) > 1.5`
#'     the data are treated as percent and divided by 100 (with a warning);
#'     otherwise passed through.
#'   \item `response_scale = "percent"`: always divide by 100.
#'   \item `response_scale = "amount"`: divide by `ll` (the larger-later
#'     reward); `ll` is required.
#' }
#' After scaling, values `> 1` are clamped to 1 and values `< 0` to 0, each
#' with a warning naming the count. `0` and `1` are valid (the SLT-beta
#' family's purpose), so they are never warned about.
#'
#' @param data A data frame containing the id, delay (x), and indifference
#'   proportion (y) columns.
#' @param y_var,x_var,id_var Column names in `data` for the response, delay,
#'   and subject id. Defaults `"y"`, `"x"`, `"id"`.
#' @param ll Optional numeric larger-later reward; required when
#'   `response_scale = "amount"`.
#' @param response_scale One of `"proportion"` (default), `"percent"`,
#'   `"amount"`.
#' @param extra_cols Optional character vector of additional column names in
#'   `data` to carry through onto the returned frame (the union of `factors`
#'   and `continuous_covariates`). These are retained verbatim so the model
#'   frame can be complete-cased once downstream; the validator never drops
#'   them. Missing names error.
#' @return A list with:
#'   \describe{
#'     \item{data}{a data frame whose first three columns are `id`, `x`, `y`,
#'       followed by any retained `extra_cols` (in the requested order).}
#'     \item{coercion_info}{list(divided_by, n_clamped_hi, n_clamped_lo,
#'       scale_detected).}
#'   }
#' @keywords internal
#' @noRd
.dd_validate_ip <- function(data,
                            y_var = "y",
                            x_var = "x",
                            id_var = "id",
                            ll = NULL,
                            extra_cols = NULL,
                            response_scale = c("proportion", "percent", "amount")) {
  response_scale <- match.arg(response_scale)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  needed <- c(id_var, x_var, y_var)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Column(s) not found in `data`: ",
      paste(shQuote(missing_cols), collapse = ", "),
      ". Set id_var / x_var / y_var to match your data.",
      call. = FALSE
    )
  }

  # Retain factor/covariate columns alongside id/x/y (row-coherence): they are
  # carried through verbatim so a single model frame can be complete-cased once
  # in .dd_tmb_prepare_data(). Drop any that collide with the canonical names.
  extra_cols <- setdiff(unique(extra_cols), c(id_var, x_var, y_var))
  if (length(extra_cols) > 0L) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra) > 0L) {
      stop(
        "Factor/covariate column(s) not found in `data`: ",
        paste(shQuote(missing_extra), collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  long <- data.frame(
    id = data[[id_var]],
    x = as.numeric(data[[x_var]]),
    y = as.numeric(data[[y_var]]),
    stringsAsFactors = FALSE
  )
  # Carry retained columns through unchanged (preserve factor class/levels).
  for (col in extra_cols) {
    long[[col]] <- data[[col]]
  }

  # --- scale detection / division ------------------------------------------
  max_y <- suppressWarnings(max(long$y, na.rm = TRUE))
  divided_by <- 1
  scale_detected <- "proportion"

  if (response_scale == "amount") {
    if (is.null(ll) || !is.finite(ll) || ll <= 0) {
      stop(
        "response_scale = 'amount' requires a positive larger-later reward `ll`.",
        call. = FALSE
      )
    }
    divided_by <- ll
    scale_detected <- "amount"
    long$y <- long$y / ll
    warning(
      sprintf(
        "Detected amount-scale responses; divided y by the larger-later reward (ll = %s).",
        format(ll)
      ),
      call. = FALSE
    )
  } else if (response_scale == "percent" ||
    (response_scale == "proportion" && is.finite(max_y) && max_y > 1.5)) {
    divided_by <- 100
    scale_detected <- "percent"
    long$y <- long$y / 100
    warning(
      "Detected percent-scale responses (max > 1.5); divided y by 100 to map to [0, 1].",
      call. = FALSE
    )
  }

  # --- clamp mild out-of-range (loud, named counts) ------------------------
  n_clamped_hi <- sum(long$y > 1, na.rm = TRUE)
  n_clamped_lo <- sum(long$y < 0, na.rm = TRUE)
  if (n_clamped_hi > 0L) long$y[!is.na(long$y) & long$y > 1] <- 1
  if (n_clamped_lo > 0L) long$y[!is.na(long$y) & long$y < 0] <- 0
  if (n_clamped_hi > 0L || n_clamped_lo > 0L) {
    warning(
      sprintf(
        "clamped out-of-range y to [0, 1]: %d value(s) > 1 set to 1, %d value(s) < 0 set to 0.",
        n_clamped_hi, n_clamped_lo
      ),
      call. = FALSE
    )
  }

  # --- error on residual NA in y -------------------------------------------
  if (anyNA(long$y)) {
    stop(
      sprintf("`%s` contains %d NA value(s) after coercion; remove or impute them.",
        y_var, sum(is.na(long$y))),
      call. = FALSE
    )
  }
  if (anyNA(long$id) || anyNA(long$x)) {
    stop("`id` and `x` must not contain NA.", call. = FALSE)
  }

  list(
    data = long,
    coercion_info = list(
      divided_by = divided_by,
      n_clamped_hi = as.integer(n_clamped_hi),
      n_clamped_lo = as.integer(n_clamped_lo),
      scale_detected = scale_detected
    )
  )
}
