#' Internal instrument registry
#'
#' Serves the per-instrument design table and ladder metadata for the
#' questionnaire scorers: the 27-item MCQ (Kirby, Petry, & Bickel, 1999),
#' the 21-item MCQ (Kirby & Marakovic, 1996), and the 30-item PDQ
#' (Madden, Petry, & Johnson, 2009).
#'
#' @param instrument One of "mcq27", "mcq21", "pdq".
#' @return List with `key`, `items`, `table` (design table in ladder order),
#'   `n_mag` (items per magnitude/block), `edge_k` (overall-ladder edge
#'   constant, or `NA` when the overall ladder uses the repeat-last
#'   convention instead), `overall_order_col` (name of the table column
#'   giving each item's position in the pooled overall ladder, or `NULL`
#'   when the table's row order already is that ladder), `value_col`,
#'   `rank_col`, `param`, and `rank_labels`.
#' @details Edge conventions: MCQ overall ladders append 0.25 (27-item) /
#'   0.1333 (21-item) past the steepest item (Kaplan et al., 2014 Excel
#'   scorers); MCQ magnitude ladders and all PDQ block ladders repeat their
#'   last indifference value (verified against all 3 x 1024 response
#'   patterns in Gray et al.'s 2016 PDQ lookup tables). The PDQ's pooled
#'   overall ladder (`overall_order_col = "overall_rank"`, repeat-last
#'   edge) is a beezdiscounting extension -- the published PDQ scoring
#'   (Madden et al., 2009; Gray et al., 2016) scores blocks only.
#' @keywords internal
.instrument_registry <- function(instrument) {
  if (
    length(instrument) != 1L ||
      !is.character(instrument) ||
      is.na(instrument) ||
      !instrument %in% c("mcq27", "mcq21", "pdq")
  ) {
    stop(
      "`instrument` must be one of \"mcq27\", \"mcq21\", \"pdq\".",
      call. = FALSE
    )
  }
  switch(
    instrument,
    mcq27 = list(
      key = "mcq27",
      items = 27L,
      table = lookup,
      n_mag = 9L,
      edge_k = 0.25,
      overall_order_col = NULL,
      value_col = "kindiff",
      rank_col = "k_rank",
      param = "k",
      rank_labels = c(
        "0.00016",
        "0.0004",
        "0.001",
        "0.0025",
        "0.006",
        "0.016",
        "0.041",
        "0.1",
        "0.25"
      )
    ),
    mcq21 = list(
      key = "mcq21",
      items = 21L,
      table = lookup21,
      n_mag = 7L,
      edge_k = 0.1333,
      overall_order_col = NULL,
      value_col = "kindiff",
      rank_col = "k_rank",
      param = "k",
      rank_labels = c(
        "0.0007",
        "0.0032",
        "0.0056",
        "0.0084",
        "0.016",
        "0.033",
        "0.13"
      )
    ),
    pdq = list(
      key = "pdq",
      items = 30L,
      table = lookup_pdq,
      n_mag = 10L,
      # no published overall-ladder edge constant: the pooled overall
      # ladder (a beezdiscounting extension) uses the repeat-last edge
      edge_k = NA_real_,
      overall_order_col = "overall_rank",
      value_col = "hindiff",
      rank_col = "h_rank",
      param = "h",
      rank_labels = as.character(1:10)
    )
  )
}

#' Internal MCQ version registry (translator)
#'
#' @param items Number of MCQ items (27 or 21)
#' @return `.instrument_registry("mcq27")` or `.instrument_registry("mcq21")`.
#' @keywords internal
.mcq_registry <- function(items) {
  if (
    length(items) != 1L ||
      !is.numeric(items) ||
      is.na(items) ||
      items != trunc(items) ||
      !items %in% c(27, 21)
  ) {
    stop("`items` must be 27 or 21.", call. = FALSE)
  }
  .instrument_registry(if (items == 27) "mcq27" else "mcq21")
}

#' Score one ascending ladder by consistency maximization
#'
#' The Kirby-style scoring core shared by the MCQ scorers (overall +
#' magnitude ladders) and the PDQ block scorer. Extracted verbatim from
#' `score_one_mcq()`'s magnitude loop; reproduces Gray et al.'s (2016)
#' published PDQ lookup tables for all 1024 response patterns per block.
#'
#' @param resp 0/1/NA responses in ladder order (1 = chose the larger
#'   delayed/risky option). Any `NA` makes all three outputs `NA` --
#'   impute upstream.
#' @param vals Ascending indifference values (k or h), one per item; must be
#'   the same length as `resp`.
#' @param edge Value appended past the steepest item (the overall-ladder
#'   edge constant, or `vals[length(vals)]` for the repeat-last convention).
#' @return `list(value, consistency, proportion)`.
#' @importFrom psych geometric.mean
#' @keywords internal
.score_ladder <- function(resp, vals, edge) {
  # the old inline code coupled these structurally; the extracted signature
  # does not, and a mismatch silently indexes past the intended ladder
  stopifnot(length(resp) == length(vals))
  n <- length(resp)
  lngth <- n + 1L
  cons <- vector(length = lngth)
  for (j in 1:lngth) {
    # first bin equals sum of larger-later/risky (1s)
    if (j == 1) {
      cons[j] <- sum(resp[j:n])
    }
    # very last bin equals sum of smaller-sooner/guaranteed (0s)
    if (j == lngth) {
      cons[j] <- sum(resp[1:j - 1] == 0)
      break()
    }
    # in between: 0s before the split plus 1s at and after it
    cons[j] <- sum(resp[1:j - 1] == 0) + sum(resp[j:n])
  }
  consmaxi <- which(cons == max(cons))
  consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
  if (0 %in% consmaxi) {
    consmaxi[which(consmaxi == 0)] <- 1
  }
  value <- if (length(consmaxi) != 0) {
    vv <- gtools::running(
      c(vals, edge)[consmaxi],
      fun = psych::geometric.mean,
      width = 2,
      by = 2
    )
    psych::geometric.mean(vv)
  } else {
    NA
  }
  list(
    value = value,
    consistency = max(cons) / n,
    proportion = cons[1] / n
  )
}
