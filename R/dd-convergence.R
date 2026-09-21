# Convergence handling shared by the TMB-backed fitters (fit_dd_tmb() and both
# fit_dd_choice() modes). Audit 2026-09-06 F-BZ4-1: a multi-start run must not
# let a non-converged start displace a converged one on NLL alone, and a
# non-converged fit must be announced at fit time and by every Wald surface,
# regardless of `verbose` (mirrors beezdemand's TICKET-063/064 gate).

#' Select the multi-start result to keep
#'
#' Candidates are ranked in three tiers: (1) converged (`opt$convergence == 0`)
#' and not flagged by `is_blowup`; (2) not flagged by `is_blowup`; (3) any
#' finite candidate. Within the first non-empty tier the lowest NLL wins. The
#' Hessian is not part of the selection (it is only computed for the winner);
#' a higher-NLL start is never substituted just to obtain standard errors.
#'
#' @param candidates List of start results (`NULL` for a failed start), each
#'   with `nll`, `opt`, and optionally `start_idx`.
#' @param is_blowup Function of `opt` returning `TRUE` for a degenerate fit.
#' @return The selected candidate with a `selection` list (`n_starts`,
#'   `n_finite`, `n_converged`, `selected_start`, `tier`,
#'   `lower_nll_nonconverged`), or `NULL` when no candidate is finite.
#' @keywords internal
#' @noRd
.dd_select_start <- function(candidates, is_blowup) {
  idx <- seq_along(candidates)
  finite <- vapply(candidates, function(r) {
    !is.null(r) && length(r$nll) == 1L && is.finite(r$nll)
  }, logical(1))
  if (!any(finite)) return(NULL)
  ok <- candidates[finite]
  ok_idx <- idx[finite]
  nll <- vapply(ok, function(r) r$nll, numeric(1))
  conv <- vapply(ok, function(r) isTRUE(r$opt$convergence == 0), logical(1))
  blow <- vapply(ok, function(r) isTRUE(is_blowup(r$opt)), logical(1))

  pick <- function(mask) {
    if (!any(mask)) return(NA_integer_)
    which(mask)[which.min(nll[mask])]
  }
  tier <- 1L
  i <- pick(conv & !blow)
  if (is.na(i)) {
    tier <- 2L
    i <- pick(!blow)
  }
  if (is.na(i)) {
    tier <- 3L
    i <- pick(rep(TRUE, length(ok)))
  }
  best <- ok[[i]]
  best$selection <- list(
    n_starts = length(candidates),
    n_finite = length(ok),
    n_converged = sum(conv),
    selected_start = best$start_idx %||% ok_idx[i],
    tier = tier,
    lower_nll_nonconverged = any(!conv & nll < nll[i])
  )
  best
}

#' Announce a non-converged fit at fit time
#' @param opt Normalized optimizer result.
#' @param what Short label of the fit ("TMB discounting fit", ...).
#' @keywords internal
#' @noRd
.dd_warn_fit_not_converged <- function(opt, what) {
  code <- opt$convergence
  msg <- opt$message %||% ""
  cli::cli_warn(
    c(
      "!" = "{what} did not converge (optimizer code {code}: {msg}).",
      "i" = "Estimates, standard errors, intervals and p-values from this fit
             may be unreliable; refit with more iterations
             ({.code tmb_control = list(iter_max = ...)}), different starts,
             or a simpler model."
    ),
    class = c("beezdiscounting_convergence_warning", "beezdiscounting_warning")
  )
  invisible(NULL)
}

#' Warn from an inference surface when the fit did not converge
#'
#' `isFALSE()` treats `NULL`/`NA` (older objects) as "no warning".
#' @param object A `beezdiscounting_tmb` or `beezdiscounting_choice` fit.
#' @keywords internal
#' @noRd
.dd_warn_if_not_converged <- function(object) {
  if (isFALSE(object$converged)) {
    cli::cli_warn(
      c(
        "!" = "This fit did not converge; its estimates, standard errors,
               intervals and p-values may be unreliable.",
        "i" = "See {.fn summary}; refit with more iterations or different
               starts before interpreting uncertainty."
      ),
      class = c("beezdiscounting_convergence_warning", "beezdiscounting_warning")
    )
  }
  invisible(NULL)
}
