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
