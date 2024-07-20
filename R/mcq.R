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
#' @return Summary dataframe
#' @export
#'
#' @examples
#' score_mcq27(mcq27)
score_mcq27 <- function(dat = dat, impute_method = "none",
                        round = 6, random = FALSE,
                        trans = "none",
                        return_data = FALSE, verbose = FALSE) {

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

    # check for 27 items
    if (!length(dat_sub$response) == 27) stop(paste0("Response length
                                                      not equal to 27
                                                      for subjectid: ", i))

    if (impute_method %in% c("inn", "INN") & any(is.na(dat_sub$response))) {
      dat_sub <- inn(dat_sub, random = random, verbose = verbose)
    }

    dfout[dfout$subjectid %in% i, 2:(ncol(dfout)-1)] <- score_one_mcq27(dat_sub,
                                                                      impute_method,
                                                                      round = round)

    if (return_data) {
      dat$newresponse[dat$subjectid == i] <- dat_sub$response
    }

  }

  dfout$impute_method <- if (!(impute_method %in% c("inn", "INN") & random)) {
    impute_method
  } else {
    "INN with random"
  }

  if (trans == "log") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_k:geomean_k, ~ log10(.x))) |>
      dplyr::rename_with(~ paste0("log10_", .x, recycle0 = TRUE),
                  overall_k:geomean_k)
  } else if (trans == "ln") {
    dfout <- dfout |>
      dplyr::mutate(dplyr::across(overall_k:geomean_k, ~ log(.x))) |>
      dplyr::rename_with(~ paste0("ln_", .x, recycle0 = TRUE),
                  overall_k:geomean_k)
  }

  if (!return_data) {
    return(dfout)
  } else {
    return(list("results" = dfout,
                "data" = dat))
  }


}

#' Score one subject's 27-item MCQ
#'
#' @param dat One subject's 27 items from the MCQ
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places
#' (passed to `base::round()`)
#'
#' @return Vector with scored 27-item MCQ metrics
#' @importFrom psych geometric.mean
#'
#' @examples
#' beezdiscounting:::score_one_mcq27(mcq27[mcq27$subjectid %in% 1, ])
score_one_mcq27 <- function(dat, impute_method = "none", round = 6) {

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

  # bring in lookup table
  dat <- merge(dat, lookup, by.x = "questionid",
               by.y = "questionid", all.x = TRUE)
  # order df
  dat <- dat[match(lookup$questionid, dat$questionid), ]

  ## overall
  # calculate consistency scores
  lngth <- 28
  cons <- vector(length = lngth)
  for (j in 1:lngth) {
    # first sorted question id equals sum of ll (1s)
    if (j == 1) cons[j] <- sum(dat$response[j:length(dat$response)])
    # very last k value bin equals sum of ss (0s)
    if (j == lngth) {
      cons[j] <- sum(dat$response[1:j - 1] == 0)
      break()
    }
    # for each question id in between, sum number of 0s before and
    # sum number of 1s at and after the current question id
    cons[j] <- sum(dat$response[1:j - 1] == 0) +
      sum(dat$response[j:length(dat$response)])
  }

  # populate consistency
  dfout["overall_consistency"] <- max(cons) / 27
  # find index where max consistency occurs
  consmaxi <- which(cons == max(cons))
  consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
  # if the highest consistency is at the first index, replace 0 index with 1
  if (0 %in% consmaxi) {
    consmaxi[which(consmaxi == 0)] <- 1
  }
  if (length(consmaxi) != 0) {
    kval <- gtools::running(c(lookup$kindiff, .25)[consmaxi],
                            fun = psych::geometric.mean, width = 2, by = 2)
    dfout["overall_k"] <- psych::geometric.mean(kval)
  } else {
    dfout["overall_k"] <- NA
  }

  dfout["overall_proportion"] <- cons[1] / 27

  for (k in 1:3) {

    dat_mag <- dat[dat$magnitude == mag[k], ]
    # calculate consistency scores
    lngth <- 10
    cons <- vector(length = lngth)
    for (j in 1:lngth) {
      # first sorted question id equals sum of ll (1s)
      if (j == 1) cons[j] <- sum(dat_mag$response[j:length(dat_mag$response)])
      # very last k value bin equals sum of ss (0s)
      if (j == lngth) {
        cons[j] <- sum(dat_mag$response[1:j - 1] == 0)
        break()
      }
      # for each question id in between, sum number of 0s before and
      # sum number of 1s at and after the current question id
      cons[j] <- sum(dat_mag$response[1:j - 1] == 0) +
        sum(dat_mag$response[j:length(dat_mag$response)])
    }

    dfout[paste0(mags[k], "_consistency")] <- max(cons) / 9
    # find index where max consistency occurs
    consmaxi <- which(cons == max(cons))
    consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
    # if the highest consistency is at the first index, replace 0 index with 1
    if (0 %in% consmaxi) {
      consmaxi[which(consmaxi == 0)] <- 1
    }
    if (length(consmaxi) != 0) {
      kval <- gtools::running(c(subset(lookup, magnitude == mag[k])$kindiff,
                                subset(lookup, magnitude == mag[k])$kindiff[10])[consmaxi],
                              fun = psych::geometric.mean, width = 2, by = 2)
      dfout[paste0(mags[k], "_k")] <- psych::geometric.mean(kval)
    } else {
      dfout[paste0(mags[k], "_k")] <- NA
    }
    dfout[paste0(mags[k], "_proportion")] <- cons[1] / 9


  }

  dfout["geomean_k"] <- psych::geometric.mean(dfout[c("small_k",
                                                      "medium_k",
                                                      "large_k")],
                              na.rm = if (impute_method %in% c("ggm", "GGM")) TRUE else FALSE)
  dfout["composite_consistency"] <- sum(dfout[c("small_consistency",
                                             "medium_consistency",
                                             "large_consistency")]) / 3

  dfout <- round(dfout, digits = round)

  return(dfout)

}

#' Calculates item nearest neighbor imputation approach discussed by
#' Yeh et al. (2023)
#'
#' @param dat A single subject's 27-item MCQ data in long form
#' @param random Boolean whether to insert a random draw (0 or 1) for NAs
#' @param verbose Boolean whether to print subject and question ids pertaining
#' to missing data
#'
#' @return An imputed data set to be scored
#'
inn <- function(dat, random, verbose) {
  dat <- merge(dat, lookup,
    by.x = "questionid",
    by.y = "questionid", all.x = TRUE
  )
  dat <- dat[match(lookup$questionid, dat$questionid), ]
  if (verbose) print(paste("NA found for id:", unique(dat$subjectid)))
  split_dfs <- suppressWarnings(split(dat, (0:nrow(dat) %/% 3))[-10])

  for (i in seq_along(split_dfs)) {
    if (!any(is.na(split_dfs[[i]]$response))) {
      next()
    } else {
      subqids <- split_dfs[[i]]$questionid
      naqids <- split_dfs[[i]]$questionid[which(is.na(split_dfs[[i]]$response))]
      if (verbose) print(paste0(c("NAs for questionids: ", naqids), collapse = " "))
      if (length(unique(split_dfs[[i]]$response[!(split_dfs[[i]]$questionid %in% naqids)])) == 1) {
        # if all non-na values are the same, replace with that non-na number
        dat$response[dat$questionid %in% naqids] <- unique(split_dfs[[i]]$response[!(split_dfs[[i]]$questionid %in% naqids)])
      } else {
        if (random) {
          dat$response[dat$questionid %in% naqids] <- sample(0:1,
            length(naqids),
            replace = TRUE
          )
        }
      }
    }
  }
  return(dat[order(as.numeric(dat$questionid)), c(2, 1, 3)])
}


#' Calculate proportion of SIR/SS responses at each k value
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#'
#' @return Dataframe with proportion of SIR/SS responses at each k rank
#' @export
#'
#' @examples prop_ss(mcq27)
prop_ss <- function(dat) {

  # bring in lookup table
  dat <- merge(dat, lookup, by.x = "questionid",
               by.y = "questionid", all.x = TRUE)
  # order df
  dat <- dat[match(lookup$questionid, dat$questionid), ]

  prop_ss_tbl <- dplyr::group_by(dat, k_rank) |>
    dplyr::summarise(prop_ss = sum(response == 0, na.rm = TRUE) / 3) |>
    dplyr::ungroup() |>
    dplyr::mutate(prop_ss = round(prop_ss, 2))

  if (any(is.na(dat$response))) {
    warning("Missing data found and ignored. Consider imputing missing data.")
  }

  return(prop_ss_tbl)

}

#' Provide a summary of the results from the MCQ ouutput table.
#'
#' @param res Dataframe with MCQ results (output from the `calc_mcq` function)
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


#' Get internal lookup table for the 27-item MCQ
#'
#' @return Dataframe with questionid, magnitude, and kindiff
#' @export
#'
#' @examples get_lookup_table()
get_lookup_table <- function() {
  return(lookup)
}
