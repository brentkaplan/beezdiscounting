# R/dd-re-utils.R

#' Normalize the `random_effects` argument to a canonical block list
#'
#' MVP scope: only a single intercept-only random effect on the discount rate,
#' supplied as the formula `k ~ 1`, is supported. The natural parameter name is
#' `"k"`; the model is fit on `log k` internally. Anything richer (random
#' slopes, intercept suppression, multiple LHS parameters, non-`k` LHS,
#' pdMat / list / pdBlocked objects) is rejected with a message pointing to the
#' fast-follow. This is the single-random-intercept special case of
#' `beezdemand`'s `.normalize_re_input()`
#' (R/random-effects-utils.R:58-115).
#'
#' @param random_effects A two-sided formula; the only supported value is
#'   `k ~ 1`.
#' @param data Optional data frame (unused in the MVP; accepted for signature
#'   compatibility with the richer fast-follow parser).
#' @return list(source = "formula", blocks = list(list(param = "k",
#'   terms = "(Intercept)", pdmat_class = "pdDiag", formula = k ~ 1, dim = 1L))).
#' @keywords internal
#' @noRd
.dd_normalize_re <- function(random_effects, data = NULL) {
  if (!inherits(random_effects, "formula")) {
    stop(
      "`random_effects` must be a formula; the only supported value in this ",
      "release is `k ~ 1` (a single random intercept on log k).",
      call. = FALSE
    )
  }
  if (length(random_effects) != 3L) {
    stop(
      "`random_effects` must be a two-sided formula `k ~ 1`.",
      call. = FALSE
    )
  }

  lhs_vars <- all.vars(random_effects[[2]])
  if (!identical(lhs_vars, "k")) {
    stop(
      "The random-effects LHS must be exactly `k` (got ",
      paste(shQuote(lhs_vars), collapse = ", "),
      "). Subject-random phi and 2-parameter random effects are out of scope ",
      "in this release.",
      call. = FALSE
    )
  }

  rhs <- random_effects[[3]]
  tt <- stats::terms(stats::as.formula(paste("~", deparse1(rhs))))
  has_intercept <- attr(tt, "intercept") == 1L
  n_terms <- length(attr(tt, "term.labels"))
  if (!has_intercept || n_terms > 0L) {
    stop(
      "Only a single intercept-only random effect `k ~ 1` is supported in ",
      "this release; random slopes are a fast-follow.",
      call. = FALSE
    )
  }

  list(
    source = "formula",
    blocks = list(list(
      param       = "k",
      terms       = "(Intercept)",
      pdmat_class = "pdDiag",
      formula     = k ~ 1,
      dim         = 1L
    ))
  )
}
