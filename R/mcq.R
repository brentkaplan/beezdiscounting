#' Score 27-item MCQ
#'
#' @param dat Dataframe (longform) with subjectid, questionid, and response
#' (0 for SIR/SS and 1 for LDR/LL)
#'
#' @return Summary dataframe
#' @export
#'
#' @examples
#' score_mcq27(mcq27)
score_mcq27 <- function(dat = dat) {

  # length of ids
  nids <- length(unique(dat$subjectid))

  # magnitudes
  mag <- c("S", "M", "L")
  mags <- c("small", "medium", "large")

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
    large_proportion = rep(NA, length = nids)
  )


  for (i in unique(dat$subjectid)) {

    # filter one subject
    dat_sub <- dat[dat$subjectid == i, ]

    # check for 27 items
    if (!length(dat_sub$response) == 27) stop(paste0("Response length not equal to 27 for subjectid: ", i))

    # bring in lookup table
    dat_sub <- merge(dat_sub, lookup, by.x = "questionid", by.y = "questionid", all.x = TRUE)
    # order df
    dat_sub <- dat_sub[match(lookup$questionid, dat_sub$questionid),]


    ## overall
    # calculate consistency scores
    lngth <- 28
    cons <- vector(length = lngth)
    for (j in 1:lngth) {
      # first sorted question id equals sum of ll (1s)
      if (j == 1) cons[j] <- sum(dat_sub$response[j:length(dat_sub$response)])
      # very last k value bin equals sum of ss (0s)
      if (j == lngth) {
        cons[j] <- sum(dat_sub$response[1:j-1] == 0)
        break()
      }
      # for each question id in between, sum number of 0s before and
      # sum number of 1s at and after the current question id
      cons[j] <- sum(dat_sub$response[1:j-1] == 0) +
        sum(dat_sub$response[j:length(dat_sub$response)])
    }

    # populate consistency
    dfout$overall_consistency[dfout$subjectid == i] <- max(cons)/27
    # find index where max consistency occurs
    consmaxi <- which(cons == max(cons))
    consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
    # if the highest consistency is at the first index, replace 0 index with 1
    if (0 %in% consmaxi) {
      consmaxi[which(consmaxi == 0)] <- 1
    }

    kval <- gtools::running(c(lookup$kindiff, .25)[consmaxi], fun = gmean, width = 2, by = 2)

    dfout$overall_k[dfout$subjectid == i] <- gmean(kval)
    dfout$overall_proportion[dfout$subjectid == i] <- cons[1]/27

    for (k in 1:3) {

      dat_sub_mag <- dat_sub[dat_sub$magnitude == mag[k], ]
      # calculate consistency scores
      lngth <- 10
      cons <- vector(length = lngth)
      for (j in 1:lngth) {
        # first sorted question id equals sum of ll (1s)
        if (j == 1) cons[j] <- sum(dat_sub_mag$response[j:length(dat_sub_mag$response)])
        # very last k value bin equals sum of ss (0s)
        if (j == lngth) {
          cons[j] <- sum(dat_sub_mag$response[1:j-1] == 0)
          break()
        }
        # for each question id in between, sum number of 0s before and
        # sum number of 1s at and after the current question id
        cons[j] <- sum(dat_sub_mag$response[1:j-1] == 0) +
          sum(dat_sub_mag$response[j:length(dat_sub_mag$response)])
      }

      dfout[[paste0(mags[k], "_consistency")]][dfout$subjectid == i] <- max(cons)/9
      # find index where max consistency occurs
      consmaxi <- which(cons == max(cons))
      consmaxi <- sort(rbind(consmaxi, (consmaxi - 1)))
      # if the highest consistency is at the first index, replace 0 index with 1
      if (0 %in% consmaxi) {
        consmaxi[which(consmaxi == 0)] <- 1
      }

      kval <- gtools::running(c(subset(lookup, magnitude == mag[k])$kindiff,
                                subset(lookup, magnitude == mag[k])$kindiff[10])[consmaxi],
                              fun = gmean, width = 2, by = 2)

      dfout[[paste0(mags[k], "_k")]][dfout$subjectid == i] <- gmean(kval)
      dfout[[paste0(mags[k], "_proportion")]][dfout$subjectid == i] <- cons[1]/9


    }

    dfout$geomean_k[dfout$subjectid == i] <- gmean(dfout[dfout$subjectid == i,
                                                         c("small_k",
                                                           "medium_k",
                                                           "large_k")])
    dfout$composite_consistency[dfout$subjectid == i] <- sum(dfout[dfout$subjectid == i,
                                                                   c("small_consistency",
                                                                     "medium_consistency",
                                                                     "large_consistency")])/3

  }

  return(dfout)

}



