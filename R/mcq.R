#' Score 27-item MCQ
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places (passed to `base::round()`)
#'
#' @return Summary dataframe
#' @export
#'
#' @examples
#' score_mcq27(mcq27)
score_mcq27 <- function(dat = dat, impute_method = "none", round = 6) {

  if (!impute_method %in% c("none", "ggm", "GGM", "inn", "INN")) {
    stop("Impute method must be one of none, ggm, GGM, inn, INN")
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
    if (!length(dat_sub$response) == 27) stop(paste0("Response length not equal to 27 for subjectid: ", i))

    dfout[dfout$subjectid %in% i, 2:(ncol(dfout)-1)] <- score_one_mcq27(dat_sub,
                                                                      impute_method,
                                                                      round = round)
  }
  dfout$impute_method <- impute_method

  return(dfout)

}


#' Score one subject's 27-item MCQ
#'
#' @param dat One subject's 27 items from the MCQ
#' @param impute_method One of: "none", "ggm", "GGM", "inn", "INN"
#' @param round Numeric specifying number of decimal places (passed to `base::round()`)
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
  # browser()
  # bring in lookup table
  dat <- merge(dat, lookup, by.x = "questionid", by.y = "questionid", all.x = TRUE)
  # order df
  dat <- dat[match(lookup$questionid, dat$questionid),]

  ## overall
  # calculate consistency scores
  lngth <- 28
  cons <- vector(length = lngth)
  for (j in 1:lngth) {
    # first sorted question id equals sum of ll (1s)
    if (j == 1) cons[j] <- sum(dat$response[j:length(dat$response)])
    # very last k value bin equals sum of ss (0s)
    if (j == lngth) {
      cons[j] <- sum(dat$response[1:j-1] == 0)
      break()
    }
    # for each question id in between, sum number of 0s before and
    # sum number of 1s at and after the current question id
    cons[j] <- sum(dat$response[1:j-1] == 0) +
      sum(dat$response[j:length(dat$response)])
  }

  # populate consistency
  dfout["overall_consistency"] <- max(cons)/27
  # find index where max consistency occurs
  consmaxi <- which(cons == max(cons))
  consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
  # if the highest consistency is at the first index, replace 0 index with 1
  if (0 %in% consmaxi) {
    consmaxi[which(consmaxi == 0)] <- 1
  }
  if (length(consmaxi) != 0) {
    kval <- gtools::running(c(lookup$kindiff, .25)[consmaxi], fun = psych::geometric.mean, width = 2, by = 2)
    dfout["overall_k"] <- psych::geometric.mean(kval)
  } else {
    dfout["overall_k"] <- NA
  }

  dfout["overall_proportion"] <- cons[1]/27

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
        cons[j] <- sum(dat_mag$response[1:j-1] == 0)
        break()
      }
      # for each question id in between, sum number of 0s before and
      # sum number of 1s at and after the current question id
      cons[j] <- sum(dat_mag$response[1:j-1] == 0) +
        sum(dat_mag$response[j:length(dat_mag$response)])
    }

    dfout[paste0(mags[k], "_consistency")] <- max(cons)/9
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
    dfout[paste0(mags[k], "_proportion")] <- cons[1]/9


  }

  dfout["geomean_k"] <- psych::geometric.mean(dfout[c("small_k", "medium_k", "large_k")],
                              na.rm = if (impute_method %in% c("ggm", "GGM")) TRUE else FALSE)
  dfout["composite_consistency"] <- sum(dfout[c("small_consistency",
                                             "medium_consistency",
                                             "large_consistency")])/3

  dfout <- round(dfout, digits = round)

  return(dfout)

}
