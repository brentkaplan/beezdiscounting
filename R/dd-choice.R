#' @useDynLib beezdiscounting, .registration = TRUE
NULL

# ==============================================================================
# Choice family (Family 2) — trial-level SS-vs-LL binomial GLMM (structural).
# ==============================================================================

#' Validate + coerce a trial-level SS-vs-LL choice data frame
#'
#' Remaps caller columns to the canonical `id, ss_amount, ll_amount, delay,
#' choice` frame, coerces types, and loudly rejects malformed data: `choice`
#' must be binary `{0,1}`; `ss_amount`/`ll_amount` finite & > 0; `delay` finite
#' & > 0 (SS is immediate, delay 0, this slice). `ll_amount <= ss_amount` warns
#' (not an error — some designs include catch trials). Declared `extra_cols`
#' (between-subject factors/covariates for the log-k design) are retained
#' verbatim (class/levels preserved).
#'
#' @param data Data frame with the choice columns.
#' @param id_var,ss_var,ll_var,delay_var,choice_var Column names.
#' @param extra_cols Optional character vector of additional columns to retain
#'   (between-subject factors/covariates).
#' @return list(data = canonical frame (+ extras), coercion_info = list()).
#' @keywords internal
#' @noRd
.dd_validate_choice <- function(data, id_var = "id", ss_var = "ss_amount",
                                ll_var = "ll_amount", delay_var = "delay",
                                choice_var = "choice", extra_cols = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  needed <- c(id_var, ss_var, ll_var, delay_var, choice_var)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(shQuote(missing_cols), collapse = ", "),
         ". Set the *_var arguments to match your data.", call. = FALSE)
  }
  long <- data.frame(
    id = as.character(data[[id_var]]),
    ss_amount = suppressWarnings(as.numeric(data[[ss_var]])),
    ll_amount = suppressWarnings(as.numeric(data[[ll_var]])),
    delay = suppressWarnings(as.numeric(data[[delay_var]])),
    choice = suppressWarnings(as.numeric(data[[choice_var]])),
    stringsAsFactors = FALSE
  )

  # choice strictly binary {0,1} (NA tolerated here; dropped in prepare).
  ch <- long$choice
  bad_choice <- !is.na(ch) & !(ch %in% c(0, 1))
  if (any(bad_choice)) {
    stop(sprintf("`%s` must be binary 0/1 (1 = chose LL); %d value(s) are not.",
                 choice_var, sum(bad_choice)), call. = FALSE)
  }

  # amounts finite & > 0.
  for (nm in c("ss_amount", "ll_amount")) {
    v <- long[[nm]]
    if (any(!is.na(v) & (!is.finite(v) | v <= 0))) {
      stop(sprintf("`%s` must be finite and positive.", nm), call. = FALSE)
    }
  }
  # delay finite & > 0.
  d <- long$delay
  if (any(!is.na(d) & !is.finite(d))) {
    stop(sprintf("`%s` must be finite.", delay_var), call. = FALSE)
  }
  if (any(!is.na(d) & d <= 0)) {
    stop(sprintf("`%s` must be > 0 (SS is immediate; LL delay positive).",
                 delay_var), call. = FALSE)
  }
  # ll > ss expected (warn).
  bad_order <- !is.na(long$ll_amount) & !is.na(long$ss_amount) &
    long$ll_amount <= long$ss_amount
  if (any(bad_order)) {
    warning(sprintf(paste0("%d row(s) have ll_amount <= ss_amount; expected the ",
                           "larger-later reward to exceed the smaller-sooner one."),
                    sum(bad_order)), call. = FALSE)
  }

  # Retain declared extra columns (between-subject design) verbatim.
  extra_cols <- setdiff(unique(extra_cols),
                        c(id_var, ss_var, ll_var, delay_var, choice_var))
  if (length(extra_cols)) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra)) {
      stop("Factor/covariate column(s) not found: ",
           paste(shQuote(missing_extra), collapse = ", "), ".", call. = FALSE)
    }
    collide <- intersect(extra_cols,
                         c("id", "ss_amount", "ll_amount", "delay", "choice"))
    if (length(collide)) {
      stop("Factor/covariate column(s) collide with reserved canonical names: ",
           paste(shQuote(collide), collapse = ", "), ".", call. = FALSE)
    }
    for (col in extra_cols) long[[col]] <- data[[col]]  # preserve class/levels
  }

  list(data = long, coercion_info = list(n_rows_in = nrow(data)))
}


#' Prepare validated choice data for TMB (0-indexed subjects, row-coherence)
#'
#' Complete-cases over all modeling columns (canonical + declared extras) ONCE,
#' builds a 0-indexed `subject_id` aligned to sorted `subject_levels` (index by
#' name, never positionally), and returns parallel arrays + the single filtered
#' frame. Retained factor extras are `droplevels()`-ed after row filtering.
#'
#' @param data Canonical frame (+ extras) from [.dd_validate_choice()].
#' @param extra_cols Optional character vector of extra columns to keep.
#' @return list(ss_amount, ll_amount, delay, choice, subject_id, subject_levels,
#'   n_subjects, n_obs, data).
#' @keywords internal
#' @noRd
.dd_choice_prepare_data <- function(data, extra_cols = NULL) {
  canonical <- c("id", "ss_amount", "ll_amount", "delay", "choice")
  model_cols <- intersect(unique(c(canonical, extra_cols)), names(data))
  keep <- stats::complete.cases(data[, model_cols, drop = FALSE])
  frame <- data[keep, , drop = FALSE]
  if (nrow(frame) == 0L) {
    stop("No complete cases remain after dropping NA rows.", call. = FALSE)
  }
  subject_levels <- sort(unique(frame$id))
  subject_map <- stats::setNames(seq_along(subject_levels) - 1L, subject_levels)
  subject_id <- as.integer(subject_map[frame$id])
  extras <- setdiff(intersect(extra_cols, names(frame)), canonical)
  cleaned <- frame[, canonical, drop = FALSE]
  for (col in extras) {
    v <- frame[[col]]
    if (is.factor(v)) v <- droplevels(v)
    cleaned[[col]] <- v
  }
  list(
    ss_amount = cleaned$ss_amount, ll_amount = cleaned$ll_amount,
    delay = cleaned$delay, choice = cleaned$choice,
    subject_id = subject_id, subject_levels = subject_levels,
    n_subjects = length(subject_levels), n_obs = nrow(cleaned),
    data = cleaned
  )
}
