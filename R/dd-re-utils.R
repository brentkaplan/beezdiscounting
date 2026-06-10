# R/dd-re-utils.R

#' Flatten an additive LHS expression into its top-level terms
#' @keywords internal
#' @noRd
.dd_re_flatten_plus <- function(e) {
  if (is.call(e) && identical(e[[1]], as.name("+"))) {
    c(.dd_re_flatten_plus(e[[2]]), .dd_re_flatten_plus(e[[3]]))
  } else {
    list(e)
  }
}

#' Parse the RE-formula LHS into canonical parameter names
#'
#' Accepts the bare symbol `k` (1-D) or the sum `k + phi` (2-D, any order).
#' Rejects transforms, non-{k,phi} symbols, a phi-only LHS, >2 params, and
#' duplicates — each pointing at the relevant fast-follow.
#' @keywords internal
#' @noRd
.dd_re_lhs_params <- function(lhs) {
  terms <- .dd_re_flatten_plus(lhs)
  syms <- vapply(terms, function(t) {
    if (!is.symbol(t)) {
      stop("The random-effects LHS must be bare symbol(s) (`k` or `k + phi`); ",
           "got `", deparse1(t), "`. Transforms of a parameter (e.g. `log(k)`) ",
           "must be a bare symbol. Subject-random `s`/`gamma`, transforms, and ",
           "random slopes are a fast-follow.", call. = FALSE)
    }
    as.character(t)
  }, character(1))
  if (anyDuplicated(syms)) {
    stop("The random-effects LHS has a duplicated parameter (got `",
         paste(syms, collapse = " + "), "`).", call. = FALSE)
  }
  bad <- setdiff(syms, c("k", "phi"))
  if (length(bad) > 0L) {
    stop("The random-effects LHS may only be `k` and/or `phi` this release ",
         "(got `", paste(bad, collapse = ", "), "`). Subject-random `s` / `gamma` ",
         "are a fast-follow.", call. = FALSE)
  }
  if (!("k" %in% syms)) {
    stop("The random-effects LHS must include `k` (a phi-only random effect is ",
         "not supported).", call. = FALSE)
  }
  intersect(c("k", "phi"), syms)   # canonical order
}

#' Normalize the `random_effects` argument to a canonical block list
#'
#' Supported: `k ~ 1` (single random intercept on log k) and `k + phi ~ 1`
#' (a correlated/independent 2-D random intercept on `(log k, log phi)`,
#' SLT-beta only — the family guard lives in `fit_dd_tmb()`). The
#' `covariance_structure` (`"pdSymm"` default, or `"pdDiag"`) sets the 2-D
#' block's class; it is irrelevant for the 1-D case. Anything richer (random
#' slopes, transformed LHS, other parameters, pdMat / pdBlocked objects) is
#' rejected with a fast-follow pointer. The canonical-block envelope mirrors
#' `beezdemand::.normalize_re_input()` so a pdMat/pdBlocked input form is a
#' clean later drop-in.
#'
#' @param random_effects A two-sided formula: `k ~ 1` or `k + phi ~ 1`.
#' @param covariance_structure `"pdSymm"` or `"pdDiag"` (2-D block only).
#' @param data Optional data frame (unused here; accepted for signature
#'   compatibility with the richer fast-follow parser).
#' @return list(source = "formula", blocks = list(list(param, terms,
#'   pdmat_class, formula, dim))). `param` is `"k"` (dim 1) or `c("k","phi")`
#'   (dim 2).
#' @keywords internal
#' @noRd
.dd_normalize_re <- function(random_effects,
                             covariance_structure = c("pdSymm", "pdDiag"),
                             data = NULL) {
  covariance_structure <- match.arg(covariance_structure)
  if (!inherits(random_effects, "formula")) {
    stop("`random_effects` must be a formula; supported values are `k ~ 1` ",
         "and `k + phi ~ 1`.", call. = FALSE)
  }
  if (length(random_effects) != 3L) {
    stop("`random_effects` must be a two-sided formula (e.g. `k ~ 1`).",
         call. = FALSE)
  }

  params <- .dd_re_lhs_params(random_effects[[2]])

  rhs <- random_effects[[3]]
  tt <- stats::terms(stats::as.formula(paste("~", deparse1(rhs))))
  has_intercept <- attr(tt, "intercept") == 1L
  n_terms <- length(attr(tt, "term.labels"))
  if (!has_intercept || n_terms > 0L) {
    stop("Only intercept-only random effects (`k ~ 1` / `k + phi ~ 1`) are ",
         "supported this release; random slopes are a fast-follow.",
         call. = FALSE)
  }

  n_params <- length(params)
  pdmat_class <- if (n_params == 1L) "pdDiag" else covariance_structure

  canon_formula <- stats::as.formula(paste(paste(params, collapse = " + "), "~ 1"))

  list(
    source = "formula",
    blocks = list(list(
      param       = params,
      terms       = "(Intercept)",
      pdmat_class = pdmat_class,
      formula     = canon_formula,
      dim         = n_params
    ))
  )
}
