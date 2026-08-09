#' Score the 30-item Probability Discounting Questionnaire (PDQ)
#'
#' Scores the Madden, Petry, & Johnson (2009) PDQ: three 10-item blocks
#' (block 1: $20 for sure vs. a chance of $80; block 2: $40 vs. $100;
#' block 3: $40 vs. $60), each an ascending ladder of h values at
#' indifference under the hyperbolic odds model `V = A / (1 + h * theta)`,
#' `theta = (1 - p) / p` (Rachlin, Raineri, & Cross, 1991). Each block is
#' scored independently by the same consistency-maximization algorithm as
#' the MCQ scorers; this implementation reproduces the Gray et al. (2016)
#' scoring syntax lookup tables for every possible response pattern.
#'
#' @param dat Dataframe (longform) with subjectid, questionid (1-30), and
#' response (0 for the smaller guaranteed reward, 1 for the larger risky
#' reward)
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#' @param random Boolean whether to insert a random draw (0 or 1) for NAs.
#' Default is FALSE
#' @param trans Transformation to apply to h values: "none", "log", or "ln".
#' Default is "none"
#' @param return_data Boolean whether to return the original data and new
#' imputed responses. Default is FALSE.
#' @param verbose Boolean whether to print subject and question ids
#' pertaining to missing data. Default is FALSE.
#'
#' @return If `return_data = FALSE` (default), a summary data frame with one
#' row per subject: pooled `overall_h` (see Details), per-block h
#' (`block1_h`, `block2_h`, `block3_h`), their arithmetic mean (`mean_h`;
#' Gray et al.'s recommended composite) and geometric mean (`geomean_h`),
#' pooled and per-block consistency plus their mean
#' (`composite_consistency`), and pooled plus per-block proportions of
#' risky choices (`block*_proportion` is Gray et al.'s risky choice ratio).
#' If `return_data = TRUE`, a list with `results` and `data` (the input plus
#' a `newresponse` column reflecting any imputation).
#' @details Each subject's data must satisfy a strict contract: exactly one
#' row per canonical question id (30 of them; no duplicates, no unknown
#' ids, none missing) and responses coded 0, 1, or `NA` (numeric, logical,
#' or character/factor values that coerce to 0/1). Malformed input errors
#' rather than silently mis-scoring. Contrast with [pdq_to_choice()]'s
#' lenient, ragged contract.
#'
#' The published scoring (Madden et al., 2009; Gray et al., 2016) has no
#' overall 30-item ladder -- blocks are scored separately, and `mean_h` is
#' Gray et al.'s recommended composite. `overall_h` and
#' `overall_consistency` are a beezdiscounting extension: all 30 items are
#' pooled into a single ascending ladder (exact-rational h order, ties
#' broken by question id; four item pairs tie exactly) and scored by the
#' same consistency-maximization algorithm with the repeat-last edge.
#'
#' INN imputation groups items sharing an h rank (one item per block).
#' Ladders with remaining `NA` responses score `NA`; Gray et al. recommend
#' excluding subjects below 80% consistency on any block.
#' @references
#' Madden, G. J., Petry, N. M., & Johnson, P. S. (2009). Pathological
#' gamblers discount probabilistic rewards less steeply than matched
#' controls. \emph{Experimental and Clinical Psychopharmacology, 17}(5),
#' 283-290. \doi{10.1037/a0016806}
#'
#' Gray, J. C., Amlung, M. T., Palmer, A. A., & MacKillop, J. (2016).
#' Syntax for calculation of discounting indices from the monetary choice
#' questionnaire and probability discounting questionnaire. \emph{Journal
#' of the Experimental Analysis of Behavior, 106}(2), 156-163.
#' \doi{10.1002/jeab.221}
#' @export
#'
#' @examples
#' score_pdq(pdq)
score_pdq <- function(
  dat = dat,
  impute_method = "none",
  round = 6,
  random = FALSE,
  trans = "none",
  return_data = FALSE,
  verbose = FALSE
) {
  reg <- .instrument_registry("pdq")

  if (!impute_method %in% c("none", "ggm", "GGM", "inn", "INN")) {
    stop("Impute method must be one of none, ggm, GGM, inn, INN")
  }

  if (!trans %in% c("none", "log", "ln")) {
    stop("Transformation must be one of 'none', 'log', 'ln'")
  }

  if (return_data) {
    dat$newresponse <- NA
  }

  nids <- length(unique(dat$subjectid))

  dfout <- data.frame(
    subjectid = unique(dat$subjectid),
    overall_h = rep(NA, length = nids),
    block1_h = rep(NA, length = nids),
    block2_h = rep(NA, length = nids),
    block3_h = rep(NA, length = nids),
    mean_h = rep(NA, length = nids),
    geomean_h = rep(NA, length = nids),
    overall_consistency = rep(NA, length = nids),
    block1_consistency = rep(NA, length = nids),
    block2_consistency = rep(NA, length = nids),
    block3_consistency = rep(NA, length = nids),
    composite_consistency = rep(NA, length = nids),
    overall_proportion = rep(NA, length = nids),
    block1_proportion = rep(NA, length = nids),
    block2_proportion = rep(NA, length = nids),
    block3_proportion = rep(NA, length = nids),
    impute_method = rep(NA, length = nids)
  )

  for (i in unique(dat$subjectid)) {
    dat_sub <- dat[dat$subjectid == i, ]

    # exact coverage of question ids 1-30: no duplicates, unknowns,
    # missing, or fractional ids (same contract as score_mcq())
    qids_raw <- suppressWarnings(as.numeric(as.character(dat_sub$questionid)))
    qids <- as.integer(qids_raw)
    non_whole_qid <- is.na(qids_raw) | qids_raw != qids
    if (
      length(qids) != reg$items ||
        anyDuplicated(qids) ||
        any(non_whole_qid) ||
        !setequal(qids[!non_whole_qid], reg$table$questionid)
    ) {
      stop(
        "Response set not equal to ",
        reg$items,
        " unique question ids for subjectid: ",
        i,
        call. = FALSE
      )
    }
    dat_sub$questionid <- qids

    # normalize responses exactly as score_mcq() does
    resp_raw <- dat_sub$response
    if (is.numeric(resp_raw)) {
      resp_num <- resp_raw
      coercion_failed <- rep(FALSE, length(resp_raw))
    } else if (is.logical(resp_raw)) {
      resp_num <- as.numeric(resp_raw)
      coercion_failed <- rep(FALSE, length(resp_raw))
    } else {
      resp_num <- suppressWarnings(as.numeric(as.character(resp_raw)))
      coercion_failed <- is.na(resp_num) & !is.na(resp_raw)
    }
    if (any(coercion_failed) || !all(resp_num %in% c(0, 1, NA))) {
      stop("Responses must be 0, 1, or NA for subjectid: ", i, call. = FALSE)
    }
    dat_sub$response <- resp_num

    if (impute_method %in% c("inn", "INN") & any(is.na(dat_sub$response))) {
      dat_sub <- inn(dat_sub, reg, random = random, verbose = verbose)
    }

    dfout[dfout$subjectid %in% i, 2:(ncol(dfout) - 1)] <- score_one_pdq(
      dat_sub,
      reg,
      impute_method,
      round = round
    )

    if (return_data) {
      idx <- which(dat$subjectid == i)
      dat$newresponse[idx] <- dat_sub$response[match(qids, dat_sub$questionid)]
    }
  }

  dfout$impute_method <- if (!(impute_method %in% c("inn", "INN") & random)) {
    impute_method
  } else {
    "INN with random"
  }

  if (trans == "log") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_h:geomean_h, ~ log10(.x))) |>
      dplyr::rename_with(
        ~ paste0("log10_", .x, recycle0 = TRUE),
        overall_h:geomean_h
      )
  } else if (trans == "ln") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_h:geomean_h, ~ log(.x))) |>
      dplyr::rename_with(
        ~ paste0("ln_", .x, recycle0 = TRUE),
        overall_h:geomean_h
      )
  }

  class(dfout) <- c("score_pdq_output", class(dfout))

  if (!return_data) {
    return(dfout)
  } else {
    return(list("results" = dfout, "data" = dat))
  }
}

#' Score one subject's PDQ
#'
#' @param dat One subject's 30 PDQ items in long form
#' @param reg Registry list from `.instrument_registry("pdq")`
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#'
#' @return Named vector with scored PDQ metrics
#' @importFrom psych geometric.mean
#' @keywords internal
score_one_pdq <- function(dat, reg, impute_method = "none", round = 6) {
  dfout <- c(
    "overall_h" = NA,
    "block1_h" = NA,
    "block2_h" = NA,
    "block3_h" = NA,
    "mean_h" = NA,
    "geomean_h" = NA,
    "overall_consistency" = NA,
    "block1_consistency" = NA,
    "block2_consistency" = NA,
    "block3_consistency" = NA,
    "composite_consistency" = NA,
    "overall_proportion" = NA,
    "block1_proportion" = NA,
    "block2_proportion" = NA,
    "block3_proportion" = NA
  )

  dat <- merge(dat, reg$table, by = "questionid", all.x = TRUE)
  dat <- dat[match(reg$table$questionid, dat$questionid), ]

  ## overall (beezdiscounting extension, not part of Madden 2009 / Gray
  ## 2016): all 30 items pooled into one ascending ladder in the
  ## precomputed pooled order (exact-rational h, questionid tie-break;
  ## reg$overall_order_col = "overall_rank"), repeat-last edge.
  ## overall_proportion (= proportion of risky choices across all 30
  ## items) falls out of the same call.
  ov_idx <- order(dat[[reg$overall_order_col]])
  vals <- dat[[reg$value_col]]
  overall <- .score_ladder(
    dat$response[ov_idx],
    vals[ov_idx],
    vals[ov_idx][reg$items]
  )
  dfout["overall_h"] <- overall$value
  dfout["overall_consistency"] <- overall$consistency
  dfout["overall_proportion"] <- overall$proportion

  for (b in 1:3) {
    dat_b <- dat[dat$block == b, ]
    hb <- reg$table[[reg$value_col]][reg$table$block == b]
    # per-block ladder, repeat-last-h edge: reproduces the Gray et al.
    # (2016) lookup tables for all 1024 patterns per block, including the
    # all-guaranteed -> maximum h and all-risky -> minimum h conventions
    res <- .score_ladder(dat_b$response, hb, hb[reg$n_mag])
    dfout[paste0("block", b, "_h")] <- res$value
    dfout[paste0("block", b, "_consistency")] <- res$consistency
    dfout[paste0("block", b, "_proportion")] <- res$proportion
  }

  na_rm <- impute_method %in% c("ggm", "GGM")
  hs <- dfout[c("block1_h", "block2_h", "block3_h")]
  dfout["mean_h"] <- mean(hs, na.rm = na_rm)
  dfout["geomean_h"] <- psych::geometric.mean(hs, na.rm = na_rm)
  dfout["composite_consistency"] <- sum(dfout[c(
    "block1_consistency",
    "block2_consistency",
    "block3_consistency"
  )]) /
    3

  dfout <- round(dfout, digits = round)

  return(dfout)
}

#' Plot PDQ Scores
#'
#' Boxplot of PDQ h metrics (block1_h, block2_h, block3_h, mean_h,
#' geomean_h, and the pooled overall_h extension), handling the same log
#' transformations as [plot.score_mcq_output()].
#'
#' @param x A data frame returned by the `score_pdq` function.
#' @param ... Additional arguments passed to methods.
#' @param xlab Label for the x-axis. Default is "Metric".
#' @param alpha Transparency of the points in the plot. Default is 0.3.
#'
#' @return A ggplot object showing the boxplot of PDQ scores.
#' @export
#'
#' @examples plot(score_pdq(pdq))
plot.score_pdq_output <- function(x, ..., xlab = "Metric", alpha = 0.3) {
  .plot_score_mcq(
    x,
    xlab = xlab,
    alpha = alpha,
    target_levels = c(
      "block1_h",
      "block2_h",
      "block3_h",
      "mean_h",
      "geomean_h",
      "overall_h"
    ),
    param = "h"
  )
}

#' Calculate proportion of guaranteed choices at each PDQ h rank
#'
#' The PDQ analog of [prop_ss()]: pools all subjects' responses and
#' reports, for each of the 10 h ranks (one item per block per rank), the
#' proportion choosing the smaller guaranteed reward. `1 - prop_sc` is the
#' risky-choice proportion at that rank.
#'
#' @param dat Dataframe (longform) with subjectid, questionid (1-30), and
#' response (0 for the guaranteed reward and 1 for the risky reward)
#'
#' @return Dataframe with proportion of guaranteed choices at each h rank
#' @export
#'
#' @examples
#' prop_sc(pdq)
prop_sc <- function(dat) {
  reg <- .instrument_registry("pdq")

  # Normalize ids the way score_pdq() does: a fractional id (1.5) or a
  # non-coercible id ("x") counts as a mismatch rather than being
  # truncated/dropped before the set comparison -- either would otherwise
  # slip past the warning while still producing an unmatched NA-rank row
  # in the merge below.
  qids_raw <- suppressWarnings(as.numeric(as.character(dat$questionid)))
  qids_int <- as.integer(qids_raw)
  bad_qid <- (is.na(qids_raw) & !is.na(dat$questionid)) |
    (!is.na(qids_raw) & qids_raw != qids_int)
  if (
    any(bad_qid) ||
      !setequal(qids_int[!is.na(qids_int) & !bad_qid], reg$table$questionid)
  ) {
    warning(
      "Observed question ids do not exactly match the 30-item PDQ design.",
      call. = FALSE
    )
  }

  dat <- merge(dat, reg$table, by = "questionid", all.x = TRUE)

  if (any(is.na(dat$response))) {
    warning("Missing data found and ignored. Consider imputing missing data.")
  }

  prop_sc_tbl <- dplyr::group_by(dat, h_rank) |>
    dplyr::summarise(prop_sc = mean(response == 0, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      prop_sc = ifelse(is.nan(prop_sc), NA_real_, round(prop_sc, 2))
    )

  class(prop_sc_tbl) <- c("prop_sc_output", class(prop_sc_tbl))
  return(prop_sc_tbl)
}

#' Plot Proportion of Guaranteed Choices by h Rank
#'
#' @param x Output from the `prop_sc` function
#' @param ... Additional arguments passed to `ggplot2::geom_point()`
#' @param pt_shape Shape of the points in the plot. Default is 21.
#' @param pt_fill Fill color of the points in the plot. Default is "white".
#' @param pt_size Size of the points in the plot. Default is 3.
#' @param title Title of the plot.
#' @param xlab Label for the x-axis. Default is "h value rank".
#' @param ylab Label for the y-axis.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples plot(prop_sc(pdq))
plot.prop_sc_output <- function(
  x,
  ...,
  pt_shape = 21,
  pt_fill = "white",
  pt_size = 3,
  title = "Proportion of guaranteed choices by h rank",
  xlab = "h value rank",
  ylab = "Proportion of guaranteed choices"
) {
  x |>
    dplyr::mutate(group = 1) |>
    ggplot2::ggplot(ggplot2::aes(
      x = factor(h_rank),
      y = prop_sc,
      group = group
    )) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      shape = pt_shape,
      fill = pt_fill,
      size = pt_size
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal()
}


#' Convert PDQ responses to a trial-level choice frame
#'
#' Reshapes long-form Probability Discounting Questionnaire (PDQ) responses
#' into a per-trial certain-versus-risky choice frame, joining each
#' `questionid` to the canonical item design (certain amount, probabilistic
#' amount, probability, and odds against `theta = (1 - p)/p`) from
#' Madden, Petry, & Johnson (2009; see [get_lookup_table()]).
#'
#' @param responses Long-form data frame with one row per PDQ item per
#'   subject, holding the columns named by `id_var`, `question_var`, and
#'   `response_var`. `response` is `0` for the smaller guaranteed reward and
#'   `1` for the larger risky reward.
#' @param id_var,question_var,response_var Column names in `responses` for
#'   the subject id, PDQ question id (1-30), and the binary choice. Defaults
#'   match the bundled `pdq` dataset (`"subjectid"`, `"questionid"`,
#'   `"response"`).
#'
#' @return A [tibble][tibble::tibble] with columns `id` (character),
#'   `sc_amount`, `lu_amount`, `prob`, `theta` (odds against winning), and
#'   `choice` (`0`/`1`, `1` = chose the risky reward), in the input row
#'   order.
#'
#' @details Unknown or non-coercible question ids raise an error rather
#'   than silently producing unmatched rows. Ragged input is allowed --
#'   subjects need not have all items -- and `NA` responses are preserved.
#'   For the strict scorer see [score_pdq()].
#'
#' @seealso [score_pdq()], [get_lookup_table()], [mcq_to_choice()]
#' @export
#'
#' @examples
#' ch <- pdq_to_choice(pdq)
#' head(ch)
pdq_to_choice <- function(
  responses,
  id_var = "subjectid",
  question_var = "questionid",
  response_var = "response"
) {
  reg <- .instrument_registry("pdq")

  if (!is.data.frame(responses)) {
    stop("`responses` must be a data frame.", call. = FALSE)
  }
  req <- c(id_var, question_var, response_var)
  missing_cols <- req[!req %in% names(responses)]
  if (length(missing_cols)) {
    stop(
      "Column(s) not found in `responses`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  id <- as.character(responses[[id_var]])

  qnum <- suppressWarnings(as.numeric(as.character(responses[[question_var]])))
  qid <- as.integer(qnum)
  i <- match(qid, reg$table$questionid)
  invalid <- is.na(i) | (!is.na(qnum) & qnum != qid)
  if (any(invalid)) {
    bad <- unique(as.character(responses[[question_var]][invalid]))
    stop(
      "Invalid questionid(s): ",
      paste(bad, collapse = ", "),
      " (must be whole numbers in 1-",
      reg$items,
      " present in the PDQ table).",
      call. = FALSE
    )
  }

  raw_choice <- responses[[response_var]]
  choice <- suppressWarnings(as.numeric(as.character(raw_choice)))
  coercion_failed <- is.na(choice) & !is.na(raw_choice)
  bad_choice <- coercion_failed | (!is.na(choice) & !(choice %in% c(0, 1)))
  if (any(bad_choice)) {
    stop(
      "`",
      response_var,
      "` must be 0/1 (1 = chose the risky reward); ",
      sum(bad_choice),
      " value(s) are not.",
      call. = FALSE
    )
  }

  tibble::tibble(
    id = id,
    sc_amount = reg$table$sc_amount[i],
    lu_amount = reg$table$lu_amount[i],
    prob = reg$table$prob[i],
    theta = reg$table$theta[i],
    choice = choice
  )
}
