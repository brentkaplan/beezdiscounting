# R/dd-validate.R

#' Validate and coerce IP-family long data to proportions in `[0, 1]`
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
#'   \item `response_scale = "proportion"` (default): the data are treated as
#'     percent and divided by 100 (with a warning) only when they are clearly
#'     percent-scaled -- a majority of the positive values exceed 1.5 and the
#'     maximum is `<= 100`. If only a few values exceed 1.5 (an ambiguous mix of
#'     proportions and out-of-range outliers) an error asks the caller to set
#'     `response_scale` explicitly, rather than silently rescaling the column.
#'     Otherwise the data pass through unchanged.
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
    # B6: reserve the canonical output names. An extra column literally named
    # id/x/y (only reachable when the role var was remapped, e.g. id_var = "pid")
    # would otherwise overwrite the canonical column in the retain loop below.
    reserved_collide <- intersect(extra_cols, c("id", "x", "y"))
    if (length(reserved_collide) > 0L) {
      stop(
        "Factor/covariate column(s) collide with the reserved canonical names ",
        "id/x/y: ", paste(shQuote(reserved_collide), collapse = ", "),
        ". Rename them (id/x/y are reserved for subject/delay/response).",
        call. = FALSE
      )
    }
  }

  # --- NA check on original y BEFORE any coercion ----------------------------
  n_na_y <- sum(is.na(data[[y_var]]))
  if (n_na_y > 0L) {
    stop(
      sprintf(
        "`%s` contains %d NA value(s); remove or impute them before validation.",
        y_var, n_na_y
      ),
      call. = FALSE
    )
  }

  long <- data.frame(
    id = data[[id_var]],
    x = suppressWarnings(as.numeric(data[[x_var]])),
    y = suppressWarnings(as.numeric(data[[y_var]])),
    stringsAsFactors = FALSE
  )
  # Carry retained columns through unchanged (preserve factor class/levels).
  for (col in extra_cols) {
    long[[col]] <- data[[col]]
  }

  # --- x must be finite and non-negative (delays cannot be Inf/NaN/negative) --
  n_bad_x <- sum(!is.finite(long$x))
  if (n_bad_x > 0L) {
    stop(
      sprintf("`%s` contains %d non-finite value(s); delays must be finite.",
        x_var, n_bad_x),
      call. = FALSE
    )
  }
  n_neg_x <- sum(long$x < 0)
  if (n_neg_x > 0L) {
    stop(
      sprintf("`%s` contains %d negative value(s); delays must be >= 0.",
        x_var, n_neg_x),
      call. = FALSE
    )
  }

  # --- y must be finite: error on Inf/-Inf/NaN (do NOT clamp Inf -> 1) -------
  # Coercion of non-numeric strings yields NA (caught by the residual-NA check
  # below); genuine Inf/-Inf/NaN values are a data error and are rejected here
  # rather than silently clamped, which would fabricate a boundary observation.
  n_inf_y <- sum(is.infinite(long$y) | is.nan(long$y))
  if (n_inf_y > 0L) {
    stop(
      sprintf("`%s` contains %d non-finite value(s) (Inf/NaN); remove them before validation.",
        y_var, n_inf_y),
      call. = FALSE
    )
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
  } else if (response_scale == "percent") {
    # Explicit request: divide by 100 as asked.
    divided_by <- 100
    scale_detected <- "percent"
    long$y <- long$y / 100
    warning(
      "Applied percent scaling as requested (response_scale = 'percent'); divided y by 100 to map to [0, 1].",
      call. = FALSE
    )
  } else if (response_scale == "proportion" && is.finite(max_y) && max_y > 1.5) {
    # B9: robust auto-detect. Divide by 100 ONLY when the data are clearly
    # percent -- a MAJORITY of the positive values exceed 1.5 AND the maximum is
    # <= 100. A few out-of-range values among otherwise-valid proportions are
    # AMBIGUOUS (e.g. a data-entry slip of 1.51 for 0.51): never silently rescale
    # a whole column from a stray value -- abort and ask for an explicit scale.
    pos <- long$y[is.finite(long$y) & long$y > 0]
    frac_big <- if (length(pos)) mean(pos > 1.5) else 0
    if (frac_big >= 0.5 && max_y <= 100) {
      divided_by <- 100
      scale_detected <- "percent"
      long$y <- long$y / 100
      warning(
        "Detected percent-scale responses (most positive values > 1.5); divided y by 100 to map to [0, 1].",
        call. = FALSE
      )
    } else {
      stop(
        sprintf(
          paste0("Ambiguous response scale: %d value(s) exceed 1.5 but the data ",
                 "are not clearly percent-scaled (max = %s). Set `response_scale` ",
                 "explicitly ('percent', 'proportion', or 'amount')."),
          sum(long$y > 1.5, na.rm = TRUE), format(max_y)
        ),
        call. = FALSE
      )
    }
  }

  # --- clamp mild out-of-range (loud, named counts) ------------------------
  n_clamped_hi <- sum(long$y > 1, na.rm = TRUE)
  n_clamped_lo <- sum(long$y < 0, na.rm = TRUE)
  if (n_clamped_hi > 0L) long$y[!is.na(long$y) & long$y > 1] <- 1
  if (n_clamped_lo > 0L) long$y[!is.na(long$y) & long$y < 0] <- 0
  if (n_clamped_hi > 0L || n_clamped_lo > 0L) {
    warning(
      sprintf(
        "Clamped out-of-range y to [0, 1]: %d value(s) > 1 set to 1, %d value(s) < 0 set to 0.",
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
