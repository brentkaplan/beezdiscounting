# R/dd-param-space.R

# NOTE: this file uses the package-level `%||%` defined once in R/utils.R.
# Do NOT redefine a local null-coalescer here.

#' Validate a report-space string
#'
#' Ported from beezdemand::beezdemand_validate_report_space
#' (R/param-space.R:20-37).
#' @keywords internal
#' @noRd
.dd_validate_report_space <- function(
  report_space,
  choices = c("natural", "log", "log10", "internal")
) {
  if (is.null(report_space) || length(report_space) != 1) {
    stop("'report_space' must be a single character value.", call. = FALSE)
  }
  report_space <- as.character(report_space)
  if (!report_space %in% choices) {
    stop(
      "'report_space' must be one of: ",
      paste(sprintf('"%s"', choices), collapse = ", "), ".",
      call. = FALSE
    )
  }
  report_space
}

#' Display label for a core term in a given report space
#'
#' Retargeted from beezdemand::beezdemand_term_display_space
#' (R/param-space.R:39-90) to the discounting core terms k / s / phi.
#' @keywords internal
#' @noRd
.dd_term_display_space <- function(term, report_space) {
  if (is.na(term) || is.null(term)) return(NA_character_)
  term <- as.character(term)
  report_space <- as.character(report_space)

  label <- function(base, suffix) {
    prefix <- switch(report_space,
      log10 = sprintf("log10(%s)", base),
      log = sprintf("log(%s)", base),
      base
    )
    paste0(prefix, suffix)
  }

  if (grepl("^k($|_|:)", term)) return(label("k", sub("^k", "", term)))
  if (grepl("^phi($|_)", term)) return(label("phi", sub("^phi", "", term)))
  if (grepl("^s($|_)", term)) return(label("s", sub("^s", "", term)))
  term
}

#' Delta-method transform of an estimate + SE between parameter spaces
#'
#' Verbatim port of beezdemand::beezdemand_transform_est_se
#' (R/param-space.R:115-178); only the function name changes.
#'
#' @param estimate,se Numeric estimate and its standard error.
#' @param from,to One of "natural", "log", "log10", "internal".
#' @return list(estimate=, se=).
#' @keywords internal
#' @noRd
.dd_transform_est_se <- function(estimate, se, from, to) {
  if (identical(from, to) || to == "internal") {
    return(list(estimate = estimate, se = se))
  }

  ln10 <- log(10)

  if (from == "natural" && to == "log10") {
    new_est <- ifelse(is.finite(estimate) & estimate > 0, log10(estimate), NA_real_)
    new_se <- ifelse(is.finite(se) & is.finite(estimate) & estimate > 0,
      se / (estimate * ln10), NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log10" && to == "natural") {
    new_est <- 10^estimate
    new_se <- ifelse(is.finite(se), ln10 * (10^estimate) * se, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log" && to == "natural") {
    new_est <- exp(estimate)
    new_se <- ifelse(is.finite(se), exp(estimate) * se, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "natural" && to == "log") {
    new_est <- ifelse(is.finite(estimate) & estimate > 0, log(estimate), NA_real_)
    new_se <- ifelse(is.finite(se) & is.finite(estimate) & estimate > 0,
      se / estimate, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log" && to == "log10") {
    new_est <- estimate / ln10
    new_se <- ifelse(is.finite(se), se / ln10, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }
  if (from == "log10" && to == "log") {
    new_est <- estimate * ln10
    new_se <- ifelse(is.finite(se), se * ln10, NA_real_)
    return(list(estimate = new_est, se = new_se))
  }

  stop("Unsupported transform from '", from, "' to '", to, "'.", call. = FALSE)
}

#' Transform the core (k/phi/s) rows of a coefficient table to a report space
#'
#' Port of beezdemand::beezdemand_transform_coef_table
#' (R/param-space.R:180-237). The core-term detector regex is retargeted to
#' `^k($|_|:)|^s($|_)|^phi($|_)` (the `:` alternative on the k branch matches
#' the renamed display terms `k:(Intercept)`, `k:condition_B`); everything else
#' (column bookkeeping, the per-row delta-method call, display labels) is unchanged.
#'
#' @keywords internal
#' @noRd
.dd_transform_coef_table <- function(
  coef_tbl,
  report_space,
  internal_space,
  term_col = "term",
  estimate_col = "estimate",
  se_col = "std.error",
  scale_col = "estimate_scale",
  display_col = "term_display"
) {
  report_space <- .dd_validate_report_space(report_space)

  if (!nrow(coef_tbl)) return(coef_tbl)
  if (!all(c(term_col, estimate_col, se_col) %in% names(coef_tbl))) {
    return(coef_tbl)
  }

  is_core <- grepl("^k($|_|:)", coef_tbl[[term_col]]) |
    grepl("^s($|_)", coef_tbl[[term_col]]) |
    grepl("^phi($|_)", coef_tbl[[term_col]])

  out <- coef_tbl
  if (!("estimate_internal" %in% names(out))) out$estimate_internal <- NA_real_
  if (!(display_col %in% names(out))) out[[display_col]] <- as.character(out[[term_col]])
  if (!(scale_col %in% names(out))) out[[scale_col]] <- NA_character_

  for (i in which(is_core)) {
    term <- as.character(out[[term_col]][i])

    from_space <- out[[scale_col]][i] %||% internal_space
    if (is.na(from_space) || !nzchar(from_space)) from_space <- internal_space

    to_space <- report_space
    if (to_space == "internal") to_space <- from_space

    trans <- .dd_transform_est_se(
      estimate = out[[estimate_col]][i],
      se = out[[se_col]][i],
      from = from_space,
      to = to_space
    )

    out$estimate_internal[i] <- out[[estimate_col]][i]
    out[[estimate_col]][i] <- trans$estimate
    out[[se_col]][i] <- trans$se

    out[[scale_col]][i] <- to_space
    out[[display_col]][i] <- .dd_term_display_space(term, to_space)
  }

  out
}
