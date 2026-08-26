#' Score an MCQ (27- or 21-item)
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#' @param items Number of MCQ items: 27 (Kirby, Petry, & Bickel, 1999) or 21
#' (Kirby & Marakovic, 1996). Default is 27.
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#' @param random Boolean whether to insert a random draw (0 or 1) for NAs.
#' Default is FALSE
#' @param return_data Boolean whether to return the original data and new
#' imputed responses. Default is FALSE.
#' @param verbose Boolean whether to print subject and question ids pertaining
#' to missing data. Default is FALSE.
#' @param trans Transformation to apply to k values: "none", "log", or "ln".
#' Default is "none"
#'
#' @return If `return_data = FALSE` (default), a summary data frame with one
#' row per subject. If `return_data = TRUE`, a list with `results` (that
#' summary data frame) and `data` (the input data plus a `newresponse`
#' column reflecting any imputation).
#' @details Each subject's data must satisfy a strict contract: exactly one
#' row per canonical question id (`items` of them; no duplicates, no unknown
#' ids, none missing) and responses coded 0, 1, or `NA` (numeric, logical, or
#' character/factor values that coerce to 0/1). Malformed input errors
#' rather than silently mis-scoring. Contrast with [mcq_to_choice()]'s
#' lenient, ragged contract, which accepts partial per-subject coverage.
#' @export
#'
#' @examples
#' score_mcq(mcq27, items = 27)
#' dat21 <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
#' score_mcq(dat21, items = 21)
score_mcq <- function(
  dat = dat,
  items = 27,
  impute_method = "none",
  round = 6,
  random = FALSE,
  trans = "none",
  return_data = FALSE,
  verbose = FALSE
) {
  reg <- .mcq_registry(items)

  if (!impute_method %in% c("none", "ggm", "GGM", "inn", "INN")) {
    stop("Impute method must be one of none, ggm, GGM, inn, INN")
  }

  if (!trans %in% c("none", "log", "ln")) {
    stop("Transformation must be one of 'none', 'log', 'ln'")
  }

  if (return_data) {
    dat$newresponse <- NA
  }

  # length of ids
  nids <- length(unique(dat$subjectid))

  # populate dataframe
  dfout <- data.frame(
    subjectid = unique(dat$subjectid),
    overall_k = rep(NA, length = nids),
    small_k = rep(NA, length = nids),
    medium_k = rep(NA, length = nids),
    large_k = rep(NA, length = nids),
    geomean_k = rep(NA, length = nids),
    overall_consistency = rep(NA, length = nids),
    small_consistency = rep(NA, length = nids),
    medium_consistency = rep(NA, length = nids),
    large_consistency = rep(NA, length = nids),
    composite_consistency = rep(NA, length = nids),
    overall_proportion = rep(NA, length = nids),
    small_proportion = rep(NA, length = nids),
    medium_proportion = rep(NA, length = nids),
    large_proportion = rep(NA, length = nids),
    impute_method = rep(NA, length = nids)
  )

  for (i in unique(dat$subjectid)) {
    # filter one subject
    dat_sub <- dat[dat$subjectid == i, ]

    # check for exact coverage of the registry's question ids: no
    # duplicates, no unknowns, none missing, and no fractional/non-whole ids
    # (as.integer(as.character()) alone would silently truncate e.g. 1.5 to
    # the valid id 1, letting a malformed id slip past coverage validation)
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

    # normalize responses the same way mcq_to_choice() does: character/factor
    # "0"/"1" score identically to numeric 0/1, and a value that fails
    # coercion (e.g. "yes") is rejected rather than silently becoming NA and
    # passing the domain check below. Numeric input is left as-is (not
    # round-tripped through as.numeric(as.character())) so its storage type
    # (e.g. integer) is preserved in `newresponse` when return_data = TRUE.
    # Logical input (TRUE/FALSE/NA) is converted explicitly via as.numeric()
    # rather than falling into the character branch below, where
    # as.numeric(as.character(TRUE)) is NA.
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

    if (impute_method %in% c("inn", "INN") && any(is.na(dat_sub$response))) {
      dat_sub <- inn(dat_sub, reg, random = random, verbose = verbose)
    }

    dfout[dfout$subjectid %in% i, 2:(ncol(dfout) - 1)] <- score_one_mcq(
      dat_sub,
      reg,
      impute_method,
      round = round
    )

    if (return_data) {
      # match on questionid within this subject's rows: dat_sub may be in a
      # different row order than dat (e.g. inn() sorts by questionid), so a
      # positional assignment here can write imputed values to the wrong
      # original-order rows when the input is not already questionid-sorted.
      idx <- which(dat$subjectid == i)
      dat$newresponse[idx] <- dat_sub$response[match(qids, dat_sub$questionid)]
    }
  }

  dfout$impute_method <- if (!(impute_method %in% c("inn", "INN") && random)) {
    impute_method
  } else {
    "INN with random"
  }

  if (trans == "log") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_k:geomean_k, ~ log10(.x))) |>
      dplyr::rename_with(
        ~ paste0("log10_", .x, recycle0 = TRUE),
        overall_k:geomean_k
      )
  } else if (trans == "ln") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_k:geomean_k, ~ log(.x))) |>
      dplyr::rename_with(
        ~ paste0("ln_", .x, recycle0 = TRUE),
        overall_k:geomean_k
      )
  }

  class(dfout) <- if (reg$items == 27L) {
    c("score_mcq27_output", class(dfout))
  } else {
    c("score_mcq_output", class(dfout))
  }

  if (!return_data) {
    return(dfout)
  } else {
    return(list("results" = dfout, "data" = dat))
  }
}

#' Score 27-item MCQ
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#' @param random Boolean whether to insert a random draw (0 or 1) for NAs.
#' Default is FALSE
#' @param return_data Boolean whether to return the original data and new
#' imputed responses. Default is FALSE.
#' @param verbose Boolean whether to print subject and question ids pertaining
#' to missing data. Default is FALSE.
#' @param trans Transformation to apply to k values: "none", "log", or "ln".
#' Default is "none"
#'
#' @return If `return_data = FALSE` (default), a summary data frame with one
#' row per subject. If `return_data = TRUE`, a list with `results` (that
#' summary data frame) and `data` (the input data plus a `newresponse`
#' column reflecting any imputation).
#' @details The subject's data must satisfy a strict contract: exactly one
#' row per canonical question id (27 of them; no duplicates, no unknown
#' ids, none missing) and responses coded 0, 1, or `NA` (numeric, logical, or
#' character/factor values that coerce to 0/1). Malformed input errors
#' rather than silently mis-scoring. Contrast with [mcq27_to_choice()]'s
#' lenient, ragged contract, which accepts partial per-subject coverage.
#' @export
#'
#' @examples
#' score_mcq27(mcq27)
score_mcq27 <- function(
  dat = dat,
  impute_method = "none",
  round = 6,
  random = FALSE,
  trans = "none",
  return_data = FALSE,
  verbose = FALSE
) {
  score_mcq(
    dat,
    items = 27,
    impute_method = impute_method,
    round = round,
    random = random,
    trans = trans,
    return_data = return_data,
    verbose = verbose
  )
}

#' Score one subject's MCQ (27- or 21-item)
#'
#' @param dat One subject's items from the MCQ
#' @param reg Registry list from `.instrument_registry()` / `.mcq_registry()`
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#'
#' @return Vector with scored MCQ metrics
#' @importFrom psych geometric.mean
#' @keywords internal
score_one_mcq <- function(dat, reg, impute_method = "none", round = 6) {
  # magnitudes
  mag <- c("S", "M", "L")
  mags <- c("small", "medium", "large")

  dfout <- c(
    "overall_k" = NA,
    "small_k" = NA,
    "medium_k" = NA,
    "large_k" = NA,
    "geomean_k" = NA,
    "overall_consistency" = NA,
    "small_consistency" = NA,
    "medium_consistency" = NA,
    "large_consistency" = NA,
    "composite_consistency" = NA,
    "overall_proportion" = NA,
    "small_proportion" = NA,
    "medium_proportion" = NA,
    "large_proportion" = NA
  )

  # bring in registry's design table
  dat <- merge(
    dat,
    reg$table,
    by = "questionid",
    all.x = TRUE
  )
  # order df
  dat <- dat[match(reg$table$questionid, dat$questionid), ]

  ## overall
  overall <- .score_ladder(dat$response, reg$table[[reg$value_col]], reg$edge_k)
  dfout["overall_k"] <- overall$value
  dfout["overall_consistency"] <- overall$consistency
  dfout["overall_proportion"] <- overall$proportion

  for (k in 1:3) {
    dat_mag <- dat[dat$magnitude == mag[k], ]
    # explicit last-k repeat == Excel scorer edge (see .score_ladder docs)
    km <- reg$table[[reg$value_col]][reg$table$magnitude == mag[k]]
    res <- .score_ladder(dat_mag$response, km, km[reg$n_mag])
    dfout[paste0(mags[k], "_k")] <- res$value
    dfout[paste0(mags[k], "_consistency")] <- res$consistency
    dfout[paste0(mags[k], "_proportion")] <- res$proportion
  }

  dfout["geomean_k"] <- psych::geometric.mean(
    dfout[c("small_k", "medium_k", "large_k")],
    na.rm = if (impute_method %in% c("ggm", "GGM")) TRUE else FALSE
  )
  dfout["composite_consistency"] <- sum(dfout[c(
    "small_consistency",
    "medium_consistency",
    "large_consistency"
  )]) /
    3

  dfout <- round(dfout, digits = round)

  return(dfout)
}

#' Score one subject's 27-item MCQ
#'
#' @param dat One subject's 27 items from the MCQ
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#'
#' @return Vector with scored 27-item MCQ metrics
#' @keywords internal
#'
#' @examples
#' beezdiscounting:::score_one_mcq27(mcq27[mcq27$subjectid %in% 1, ])
score_one_mcq27 <- function(dat, impute_method = "none", round = 6) {
  score_one_mcq(dat, .mcq_registry(27), impute_method, round)
}

#' Calculates item nearest neighbor imputation approach discussed by
#' Yeh et al. (2023)
#'
#' @param dat A single subject's MCQ data in long form
#' @param reg Registry list from `.instrument_registry()` / `.mcq_registry()`
#' @param random Boolean whether to insert a random draw (0 or 1) for NAs
#' @param verbose Boolean whether to print subject and question ids pertaining
#' to missing data
#'
#' @return An imputed data set to be scored
#'
#' @keywords internal
inn <- function(dat, reg, random, verbose) {
  dat <- merge(dat, reg$table, by = "questionid", all.x = TRUE)
  dat <- dat[match(reg$table$questionid, dat$questionid), ]
  if (verbose) {
    print(paste("NA found for id:", unique(dat$subjectid)))
  }
  ## Yeh et al. (2023) nearest neighbors = items sharing a rank on the
  ## instrument's ladder. For the 27-item MCQ these are the historical
  ## consecutive k_rank triples; for the 21-item MCQ group sizes vary
  ## (3,3,2,4,3,3,3); for the PDQ each h_rank group holds one item per
  ## block. Upstream validation in score_mcq() guarantees canonical
  ## questionid coverage, so no NA ranks.
  split_dfs <- split(
    dat,
    factor(dat[[reg$rank_col]], levels = unique(dat[[reg$rank_col]]))
  )

  for (i in seq_along(split_dfs)) {
    if (!any(is.na(split_dfs[[i]]$response))) {
      next()
    } else {
      naqids <- split_dfs[[i]]$questionid[which(is.na(split_dfs[[i]]$response))]
      if (verbose) {
        print(paste0(c("NAs for questionids: ", naqids), collapse = " "))
      }
      if (
        length(unique(split_dfs[[i]]$response[
          !(split_dfs[[i]]$questionid %in% naqids)
        ])) ==
          1
      ) {
        # if all non-na values are the same, replace with that non-na number
        dat$response[dat$questionid %in% naqids] <- unique(split_dfs[[
          i
        ]]$response[!(split_dfs[[i]]$questionid %in% naqids)])
      } else {
        if (random) {
          dat$response[dat$questionid %in% naqids] <- sample(
            0:1,
            length(naqids),
            replace = TRUE
          )
        }
      }
    }
  }
  dat <- dat[order(as.numeric(dat$questionid)), ]
  dat[, c("subjectid", "questionid", "response")]
}


#' Calculate proportion of SIR/SS responses at each k value
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#' @param items Number of MCQ items (27 or 21)
#'
#' @details `items` must match the instrument actually administered.
#' Question ids 1-21 are valid in both the 21- and 27-item designs, so
#' passing the wrong `items` does not error; it silently pools responses
#' into the wrong k-rank rows. If the observed question ids do not exactly
#' match the requested design, `prop_ss()` warns.
#'
#' @return Dataframe with proportion of SIR/SS responses at each k rank
#' @export
#'
#' @examples
#' prop_ss(mcq27)
#' dat21 <- data.frame(subjectid = 1, questionid = 1:21, response = 1)
#' prop_ss(dat21, items = 21)
prop_ss <- function(dat, items = 27) {
  reg <- .mcq_registry(items)

  # warn (do not error -- prop_ss() has no ragged contract to protect) when
  # the observed question ids don't exactly cover this design's canonical
  # set: question ids 1-21 are valid in both the 21- and 27-item designs, so
  # e.g. prop_ss(mcq21) with the default items = 27 would otherwise silently
  # pool responses into the wrong rank table.
  observed_qids <- suppressWarnings(as.integer(as.character(dat$questionid)))
  if (!setequal(observed_qids[!is.na(observed_qids)], reg$table$questionid)) {
    warning(
      "Observed question ids do not exactly match the ", reg$items,
      "-item MCQ design (items = ", reg$items, "). `items` must match the ",
      "instrument actually administered -- a mismatch silently pools ",
      "responses under the wrong k ranks.",
      call. = FALSE
    )
  }

  # bring in lookup table (k_rank etc.); keep every row so all respondents pool
  dat <- merge(dat, reg$table, by.x = "questionid",
               by.y = "questionid", all.x = TRUE)

  if (any(is.na(dat$response))) {
    warning("Missing data found and ignored. Consider imputing missing data.")
  }

  prop_ss_tbl <- dplyr::group_by(dat, k_rank) |>
    dplyr::summarise(prop_ss = mean(response == 0, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(prop_ss = ifelse(is.nan(prop_ss), NA_real_, round(prop_ss, 2)))

  if (reg$items != 27L) {
    attr(prop_ss_tbl, "mcq_items") <- 21L
  }

  class(prop_ss_tbl) <- c("prop_ss_output", class(prop_ss_tbl))
  return(prop_ss_tbl)

}


#' Provide a summary of the results from the MCQ output table.
#'
#' @param res Dataframe with MCQ results (output from `score_mcq()` or
#' `score_mcq27()`)
#' @param na.rm Boolean whether to remove NAs from the calculation
#'
#' @return Dataframe with summary statistics
#' @export
#'
#' @examples summarize_mcq(score_mcq27(mcq27))
summarize_mcq <- function(res, na.rm = TRUE) {
  sum_tab <- res %>%
    dplyr::summarise(
      dplyr::across(dplyr::contains("overall_k"):composite_consistency, list(
        Mean = ~ mean(., na.rm = na.rm),
        SD = ~ sd(., na.rm = na.rm),
        SEM = ~ sd(., na.rm = na.rm) / sqrt(dplyr::n())
      ), .names = "{.col}-{.fn}")
    ) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = c(".value", "Statistic"), names_sep = "-") %>%
    tidyr::pivot_longer(-Statistic, names_to = "Metric", values_to = "value") |>
    tidyr::pivot_wider(names_from = Statistic, values_from = value)

  return(sum_tab)
}


#' Get an internal instrument design table
#'
#' @param items **Deprecated.** Number of MCQ items (27 or 21), kept as a
#'   back-compatible alias for `instrument = "mcq27"` / `"mcq21"`
#'   (including positional calls like `get_lookup_table(21)`). Use
#'   `instrument` in new code; `items` cannot address the PDQ.
#' @param instrument Instrument key: `"mcq27"` (Kirby, Petry, & Bickel,
#'   1999), `"mcq21"` (Kirby & Marakovic, 1996), or `"pdq"` (Madden,
#'   Petry, & Johnson, 2009). Supply `items` or `instrument`, not both.
#'   The default (neither) returns the 27-item MCQ design.
#' @return Data frame with the complete item design. MCQ tables:
#'   questionid, magnitude, kindiff, k_rank, ss_amount, ll_amount, delay
#'   (days). PDQ table: questionid, block, h_rank, overall_rank (the
#'   item's position in the pooled overall ladder used by the
#'   `overall_h` extension in [score_pdq()]), sc_amount, lu_amount,
#'   prob, theta (odds against), hindiff.
#' @export
#'
#' @examples
#' get_lookup_table()
#' get_lookup_table(items = 21)
#' get_lookup_table(instrument = "pdq")
get_lookup_table <- function(items = NULL, instrument = NULL) {
  if (!is.null(items) && !is.null(instrument)) {
    stop("Supply `items` or `instrument`, not both.", call. = FALSE)
  }
  if (!is.null(instrument)) {
    return(.instrument_registry(instrument)$table)
  }
  .mcq_registry(items %||% 27)$table
}


#' Convert 27- or 21-item MCQ responses to a trial-level choice frame
#'
#' Reshapes long-form Monetary Choice Questionnaire (MCQ) responses into
#' the per-trial smaller-sooner versus larger-later choice frame consumed by
#' [fit_dd_choice()], joining each `questionid` to the canonical item design
#' (immediate amount, delayed amount, delay) bundled in the lookup table (see
#' [get_lookup_table()]; Kirby, Petry, & Bickel (1999) for the 27-item MCQ,
#' or Kirby & Maraković (1996, Table 1) for the 21-item MCQ).
#'
#' @param responses Long-form data frame with one row per MCQ item per subject,
#'   holding the columns named by `id_var`, `question_var`, and `response_var`.
#'   `response` is `0` for the smaller-immediate reward (SIR/SS) and `1` for the
#'   larger-delayed reward (LDR/LL), the same coding [fit_dd_choice()] expects,
#'   so no recoding is applied.
#' @param items Number of MCQ items: 27 (Kirby, Petry, & Bickel, 1999) or 21
#'   (Kirby & Maraković, 1996). Default is 27.
#' @param id_var,question_var,response_var Column names in `responses` for the
#'   subject id, MCQ question id, and the binary choice. Defaults match the
#'   bundled `mcq27` dataset (`"subjectid"`, `"questionid"`, `"response"`).
#'
#' @return A [tibble][tibble::tibble] with columns `id` (character), `ss_amount`,
#'   `ll_amount`, `delay` (days), and `choice` (`0`/`1`, `1` = chose LL), in the
#'   input row order. Ready to pass to [fit_dd_choice()].
#'
#' @details Unknown or non-coercible question ids raise an error rather than
#'   silently producing unmatched rows. Question ids 1-21 are valid in both
#'   the 21- and 27-item designs, so passing the wrong `items` does not
#'   error; it silently returns the wrong amounts/delays for those ids.
#'   Make sure `items` matches the instrument actually administered. Ragged
#'   input is allowed (subjects need not have all items) and `NA`
#'   responses are preserved (they are complete-cased by [fit_dd_choice()]).
#'   For the strict scorer see [score_mcq()].
#'
#' @seealso [fit_dd_choice()], [score_mcq()], [get_lookup_table()],
#'   [mcq27_to_choice()]
#' @export
#'
#' @examples
#' ch <- mcq_to_choice(mcq27, items = 27)
#' head(ch)
#' # feeds directly into the structural choice model (requires TMB):
#' # fit_dd_choice(ch, mode = "structural", equation = "mazur")
mcq_to_choice <- function(responses,
                          items = 27,
                          id_var = "subjectid",
                          question_var = "questionid",
                          response_var = "response") {
  reg <- .mcq_registry(items)

  if (!is.data.frame(responses)) {
    stop("`responses` must be a data frame.", call. = FALSE)
  }
  req <- c(id_var, question_var, response_var)
  missing_cols <- req[!req %in% names(responses)]
  if (length(missing_cols)) {
    stop("Column(s) not found in `responses`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  id <- as.character(responses[[id_var]])

  # Coerce the question id to a whole-number key (character/factor/numeric all
  # normalize). A non-coercible, out-of-range, or fractional value (e.g. 1.5,
  # which as.integer() would silently truncate to 1) is rejected, not mapped.
  qnum <- suppressWarnings(as.numeric(as.character(responses[[question_var]])))
  qid <- as.integer(qnum)
  i <- match(qid, reg$table$questionid)
  invalid <- is.na(i) | (!is.na(qnum) & qnum != qid)
  if (any(invalid)) {
    bad <- unique(as.character(responses[[question_var]][invalid]))
    stop("Invalid questionid(s): ", paste(bad, collapse = ", "),
         " (must be whole numbers in 1-", reg$items,
         " present in the MCQ table).",
         call. = FALSE)
  }

  # `response` already uses the choice coding (1 = LL); coerce via character so
  # factor levels "0"/"1" map to the labels, not the underlying integer codes.
  raw_choice <- responses[[response_var]]
  choice <- suppressWarnings(as.numeric(as.character(raw_choice)))
  # Flag values outside {0, 1}, or that became NA only through a failed coercion
  # (e.g. "x"); a genuine NA in the input is preserved (complete-cased downstream).
  coercion_failed <- is.na(choice) & !is.na(raw_choice)
  bad_choice <- coercion_failed | (!is.na(choice) & !(choice %in% c(0, 1)))
  if (any(bad_choice)) {
    stop("`", response_var, "` must be 0/1 (1 = chose LL); ",
         sum(bad_choice), " value(s) are not.", call. = FALSE)
  }

  tibble::tibble(
    id = id,
    ss_amount = reg$table$ss_amount[i],
    ll_amount = reg$table$ll_amount[i],
    delay = reg$table$delay[i],
    choice = choice
  )
}


#' Convert 27-item MCQ responses to a trial-level choice frame
#'
#' @inherit mcq_to_choice params return details description
#' @param id_var,question_var,response_var Column names in `responses` for the
#'   subject id, MCQ question id (1-27), and the binary choice. Defaults match the
#'   bundled `mcq27` dataset (`"subjectid"`, `"questionid"`, `"response"`).
#'
#' @seealso [fit_dd_choice()], [score_mcq27()], [get_lookup_table()],
#'   [mcq_to_choice()]
#' @export
#'
#' @examples
#' ch <- mcq27_to_choice(mcq27)
#' head(ch)
#' # feeds directly into the structural choice model (requires TMB):
#' # fit_dd_choice(ch, mode = "structural", equation = "mazur")
mcq27_to_choice <- function(responses,
                            id_var = "subjectid",
                            question_var = "questionid",
                            response_var = "response") {
  mcq_to_choice(responses, items = 27, id_var = id_var,
                question_var = question_var, response_var = response_var)
}


#' Plot Proportion of SIR/SS Choices by k Value
#'
#' Plots the proportion of SIR/SS
#' choices by k value using the output of the `prop_ss` function.
#'
#' @param x Output from the `prop_ss` function
#' @param ... Additional arguments passed to `ggplot2::geom_point()`
#' @param pt_shape Shape of the points in the plot. Default is 21.
#' @param pt_fill Fill color of the points in the plot. Default is "white".
#' @param pt_size Size of the points in the plot. Default is 3.
#' @param title Title of the plot. Default is "Proportion of SIR/SS choices by k value".
#' @param xlab Label for the x-axis. Default is "k value rank".
#' @param ylab Label for the y-axis. Default is "Proportion of SS choices".
#'
#' @return A ggplot object.
#' @export
#'
#' @examples plot(prop_ss(mcq27))
plot.prop_ss_output <- function(
    x,
    ...,
    pt_shape = 21,
    pt_fill = "white",
    pt_size = 3,
    title = "Proportion of SIR/SS choices by k value",
    xlab = "k value rank",
    ylab = "Proportion of SS choices"
    ) {
  # index (not append) the labels by the k_rank values actually present, so
  # a missing rank (first or middle) can't shift labels onto the wrong
  # remaining ranks -- factor(k_rank)'s levels are these same sorted values,
  # in the same order.
  levels_present <- sort(unique(x$k_rank))
  labs_k <- .mcq_registry(attr(x, "mcq_items") %||% 27L)$rank_labels[levels_present]
  x |>
    dplyr::mutate(group = 1) |>
    ggplot2::ggplot(ggplot2::aes(x = factor(k_rank), y = prop_ss, group = group)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      ...,
      shape = pt_shape,
      fill = pt_fill,
      size = pt_size
    ) +
    ggplot2::labs(title = title,
         x = xlab,
         y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(labels = labs_k)
}


#' Plot questionnaire scores (internal implementation)
#'
#' Shared boxplot implementation for `plot.score_mcq27_output()`,
#' `plot.score_mcq_output()`, and `plot.score_pdq_output()`.
#'
#' @param x A data frame returned by `score_mcq()`, `score_mcq27()`, or
#'   `score_pdq()`.
#' @param xlab Label for the x-axis.
#' @param alpha Transparency of the points in the plot.
#' @param target_levels Metric column names, in display order.
#' @param param Parameter letter used in column suffixes and the y label
#'   ("k" or "h").
#'
#' @return A ggplot object showing the boxplot of scores.
#' @keywords internal
.plot_score_mcq <- function(
  x,
  xlab = "Metric",
  alpha = 0.3,
  target_levels = c("small_k", "medium_k", "large_k", "geomean_k", "overall_k"),
  param = "k"
) {
  suffix <- paste0("_", param)
  tmp <- x |>
    dplyr::select(dplyr::contains(c("id", suffix))) |>
    tidyr::pivot_longer(
      cols = dplyr::contains(suffix),
      names_to = "metric",
      values_to = "value"
    )

  if (any(grepl("log10", tmp$metric))) {
    tmp <- tmp |>
      dplyr::mutate(
        metric = gsub("log10_", "", metric),
        metric = factor(metric, levels = target_levels)
      )
    ylab = paste0("Log10(", param, ") value")
  } else if (any(grepl("ln", tmp$metric))) {
    tmp <- tmp |>
      dplyr::mutate(
        metric = gsub("ln_", "", metric),
        metric = factor(metric, levels = target_levels)
      )
    ylab <- paste0("Ln(", param, ") value")
  } else {
    tmp <- tmp |>
      dplyr::mutate(metric = factor(metric, levels = target_levels))
    ylab <- paste0(param, " value")
  }

  bplot <- tmp %>%
    ggplot2::ggplot(ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = xlab,
      y = ylab
    )

  return(bplot)
}

#' Plot MCQ-27 Scores
#'
#' Boxplots the MCQ-27 scores for
#' different metrics (small_k, medium_k, large_k, geomean_k, overall_k).
#' Log-transformed k values (`trans = "log"`/`"ln"`) are detected
#' and the y-axis label adjusted accordingly.
#'
#' @param x A data frame returned by the `score_mcq27` function.
#' @param ... Additional arguments passed to methods.
#' @param xlab Label for the x-axis. Default is "Metric".
#' @param alpha Transparency of the points in the plot. Default is 0.3.
#'
#' @return A ggplot object showing the boxplot of MCQ-27 scores.
#' @export
#'
#' @examples plot(score_mcq27(mcq27))
plot.score_mcq27_output <- function(x, ..., xlab = "Metric", alpha = 0.3) {
  .plot_score_mcq(x, xlab = xlab, alpha = alpha)
}

#' Plot MCQ Scores
#'
#' Boxplots the MCQ scores for
#' different metrics (small_k, medium_k, large_k, geomean_k, overall_k).
#' Log-transformed k values (`trans = "log"`/`"ln"`) are detected
#' and the y-axis label adjusted accordingly.
#'
#' @param x A data frame returned by the `score_mcq` function.
#' @param ... Additional arguments passed to methods.
#' @param xlab Label for the x-axis. Default is "Metric".
#' @param alpha Transparency of the points in the plot. Default is 0.3.
#'
#' @return A ggplot object showing the boxplot of MCQ scores.
#' @export
#'
#' @examples
#' dat21 <- data.frame(subjectid = 1:5, questionid = rep(1:21, 5), response = 1)
#' plot(score_mcq(dat21, items = 21))
plot.score_mcq_output <- function(x, ..., xlab = "Metric", alpha = 0.3) {
  .plot_score_mcq(x, xlab = xlab, alpha = alpha)
}
