#' Score 5.5 trial delay discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with id, indexes, response, k value, and effective delay 50.
#' @details
#' Currently assumes the attending questions are present and labeled "Attend-LL" and "Attend-SS"
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' score_dd(five.fivetrial_dd)
score_dd <- function(df) {
  dd1 <- df |>
    dplyr::select(ResponseId, paste0("I", seq(1, 31, by = 2)), "AttendSS" = `Attend-SS`,
                  "AttendLL" = `Attend-LL`) |>
    dplyr::select(-dplyr::contains("Timing"), -dplyr::contains("DO"))
  ddframe <- dd1 |>
    tidyr::pivot_longer(cols = 2:ncol(dd1), names_to = "index", values_to = "response") %>%
    dplyr::filter(complete.cases(.)) |>
    dplyr::mutate(response = ifelse(stringr::str_detect(response, "now"), "ss", "ll"))
  ddframe$kval <- NA
  ddframe$attentionflag <- "No"
  indexes <- paste0("I", seq(1, 31, by = 2))
  recodess <- c("24", "9.797958971", "4.898979486", "2.309401077", "0.816496581", "0.40824829",
                "0.188982237", "0.08247861", "0.039551962", "0.013411573", "0.006705787",
                "0.003353031", "0.001117723", "0.000612202", "0.000279431", "0.000129064")
  recodell <- c("16.97056275", "6.92820323", "3.265986324", "1.414213562", "0.577350269", "0.288675135",
                "0.116642369", "0.058321184", "0.023229526", "0.009483414", "0.004741707",
                "0.001935953", "0.000790349", "0.000432892", "0.000186287", "0.000109514")
  for (i in seq_along(indexes)) {
    if (length(ddframe$response[ddframe$index == indexes[i]]) == 0) next
    ddframe$kval[ddframe$index == indexes[i]] <- ifelse (ddframe$response[ddframe$index == indexes[i]] %in% "ss",
                                                            recodess[i], recodell[i])
  }
  subset_indices <- which(ddframe$ResponseId %in% ddframe$ResponseId[ddframe$index %in% "AttendSS"] & ddframe$index %in% "I1")
  ddframe <- if (length(subset_indices) != 0) ddframe[-subset_indices, ] else ddframe
  subset_indices <- which(ddframe$ResponseId %in% ddframe$ResponseId[ddframe$index %in% "AttendLL"] & ddframe$index %in% "I31")
  ddframe <- if (length(subset_indices) != 0) ddframe[-subset_indices, ] else ddframe
  ddframe$attentionflag[ddframe$index %in% "AttendSS" & ddframe$response %in% "ss"] <- "Yes"
  ddframe$attentionflag[ddframe$index %in% "AttendLL" & ddframe$response %in% "ll"] <- "Yes"
  ddframe$attentionflag[ddframe$ResponseId %in% ddframe$ResponseId[which(ddframe$attentionflag == "Yes")]]  <- "Yes"
  ddframe$kval[ddframe$index %in% "AttendSS" & ddframe$response %in% "ll"] <- "24"
  ddframe$kval[ddframe$index %in% "AttendLL" & ddframe$response %in% "ss"] <- "0.000109514"

  ddframe$kval <- as.numeric(ddframe$kval)
  ddframe$ed50 <- 1/ddframe$kval
  return(ddframe)

}

#' Extract timing metrics from 5.5 trial delay discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with ResponseId, indexes, values and timing
#' @details
#' Currently assumes the attending questions are present and labeled "Attend-LL" and "Attend-SS"
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' timing_dd(five.fivetrial_dd)
timing_dd <- function(df) {
  timing <- dplyr::select(df, ResponseId, dplyr::contains("Timing"))
  colnames(timing) <- gsub("Timing_First Click", "firstclick", colnames(timing))
  colnames(timing) <- gsub("Timing_Last Click", "lastclick", colnames(timing))
  colnames(timing) <- gsub("Timing_Page Submit", "pagesubmit", colnames(timing))
  colnames(timing) <- gsub("Timing_Click Count", "totalclicks", colnames(timing))
  timing <- timing |>
    tidyr::pivot_longer(cols = 2:ncol(timing), names_to = "question", values_to = "value") %>%
    dplyr::filter(complete.cases(.))
  timing$q <- NA
  timing$question <- gsub("Attend-LL", "AttendLL", timing$question)
  timing$question <- gsub("Attend-SS", "AttendSS", timing$question)
  timing <- timing |>
    tidyr::separate(question, c("index", "measure"), sep = "-") |>
    tidyr::spread(measure, value)
  timing$q[timing$index %in% c("I16")] <- 1
  timing$q[timing$index %in% c("I8", "I24")] <- 2
  timing$q[timing$index %in% c("I4", "I12", "I20", "I28")] <- 3
  timing$q[timing$index %in% c("I2", "I6", "I10", "I14", "I18", "I22", "I26", "I30")] <- 4
  timing$q[timing$index %in% c("I1", "I3", "I5", "I7", "I9", "I11", "I13", "I15", "I17", "I19",
                         "I21", "I23", "I25", "I27", "I29", "I31")] <- 5
  timing$q[timing$index %in% c("AttendSS", "AttendLL")] <- 6
  timing[4:7] <- sapply(timing[4:7], as.numeric)
  return(dplyr::arrange(timing, ResponseId, q))

}

#' Converts answers from 5.5 trial delay discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with the ResponseId, index, and response (ss or ll).
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' ans_dd(five.fivetrial_dd)
ans_dd <- function(df) {
  ans <- df |>
    dplyr::select(ResponseId, paste0("I", 1:31), dplyr::starts_with("Attend")) |>
    dplyr::select(-dplyr::contains("Timing"), -dplyr::contains("_DO")) %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "index", values_to = "response") %>%
    dplyr::filter(complete.cases(.)) |>
    dplyr::mutate(response = ifelse(stringr::str_detect(response, "now"), "ss", "ll"))
  ans$index <- gsub("-", "", ans$index)
  return(ans)
}

#' Calculate scores, answers, and timing for 5.5 trial delay discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns from the template.
#'
#' @return A dataframe with k/ed50 values, answers, timing
#' @export
#'
#' @examples
#' calc_dd(five.fivetrial_dd)
calc_dd <- function(df) {
    return(dplyr::left_join(timing_dd(df), ans_dd(df), by = c("ResponseId", "index")) |>
           dplyr::left_join(dplyr::select(score_dd(df), ResponseId, attentionflag, kval, ed50),
                            by = c("ResponseId")) |>
           dplyr::arrange(ResponseId, q))
}

#' Score 5.5 trial probability discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with id, indexes, response, h value, and effective probability 50.
#' @details
#' Currently assumes the attending questions are present and labeled "Attend-LL" and "Attend-SS"
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' score_pd(five.fivetrial_pd)
score_pd <- function(df) {
  pd1 <- df |>
    dplyr::select(ResponseId, paste0("I", seq(1, 31, by = 2)), "AttendSS" = `Attend-SS`,
                  "AttendLL" = `Attend-LL`) |>
    dplyr::select(-dplyr::contains("Timing"), -dplyr::contains("DO"))
  pdframe <- pd1 |>
    tidyr::pivot_longer(cols = 2:ncol(pd1), names_to = "index", values_to = "response") %>%
    dplyr::filter(complete.cases(.)) |>
    dplyr::mutate(response = ifelse(stringr::str_detect(response, "for sure"), "sc", "lu"))
  pdframe$hval <- NA
  pdframe$attentionflag <- "No"
  indexes <- paste0("I", seq(1, 31, by = 2))
  recodesc <- c("99", "56.72448", "35.510562", "19.390719", "10.173495", "5.290003",
                "2.934058", "1.601445", "0.850963", "0.457604", "0.256064",
                "0.134491", "0.074501", "0.035898", "0.022875", "0.012403")
  recodelu <- c("80.628779", "43.714986", "27.856777", "13.422618", "7.435436", "3.905279",
                "2.185294", "1.175139", "0.624436", "0.340825", "0.189036",
                "0.098295", "0.051571", "0.028161", "0.017629", "0.010101")
  for (i in seq_along(indexes)) {
    if (length(pdframe$response[pdframe$index == indexes[i]]) == 0) next
    pdframe$hval[pdframe$index == indexes[i]] <- ifelse (pdframe$response[pdframe$index == indexes[i]] %in% "sc",
                                                         recodesc[i], recodelu[i])
  }
  subset_indices <- which(pdframe$ResponseId %in% pdframe$ResponseId[pdframe$index %in% "AttendSS"] & pdframe$index %in% "I1")
  pdframe <- if (length(subset_indices) != 0) pdframe[-subset_indices, ] else pdframe
  subset_indices <- which(pdframe$ResponseId %in% pdframe$ResponseId[pdframe$index %in% "AttendLL"] & pdframe$index %in% "I31")
  pdframe <- if (length(subset_indices) != 0) pdframe[-subset_indices, ] else pdframe
  pdframe$attentionflag[pdframe$index %in% "AttendSS" & pdframe$response %in% "sc"] <- "Yes"
  pdframe$attentionflag[pdframe$index %in% "AttendLL" & pdframe$response %in% "lu"] <- "Yes"
  pdframe$attentionflag[pdframe$ResponseId %in% pdframe$ResponseId[which(pdframe$attentionflag == "Yes")]]  <- "Yes"
  pdframe$hval[pdframe$index %in% "AttendSS" & pdframe$response %in% "lu"] <- "99"
  pdframe$hval[pdframe$index %in% "AttendLL" & pdframe$response %in% "sc"] <- "0.010101"

  pdframe$hval <- as.numeric(pdframe$hval)
  pdframe$etheta50 <- 1/pdframe$hval
  pdframe$ep50 <- (1 / (pdframe$etheta50 + 1)) * 100
  return(pdframe)

}

#' Extract timing metrics from 5.5 trial probability discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with ResponseId, indexes, values and timing
#' @details
#' Currently assumes the attending questions are present and labeled "Attend-LL" and "Attend-SS"
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' timing_pd(five.fivetrial_pd)
timing_pd <- function(df) {
  timing <- dplyr::select(df, ResponseId, dplyr::contains("Timing"))
  colnames(timing) <- gsub("Timing_First Click", "firstclick", colnames(timing))
  colnames(timing) <- gsub("Timing_Last Click", "lastclick", colnames(timing))
  colnames(timing) <- gsub("Timing_Page Submit", "pagesubmit", colnames(timing))
  colnames(timing) <- gsub("Timing_Click Count", "totalclicks", colnames(timing))
  timing <- timing |>
    tidyr::pivot_longer(cols = 2:ncol(timing), names_to = "question", values_to = "value") %>%
    dplyr::filter(complete.cases(.))
  timing$q <- NA
  timing$question <- gsub("Attend-LL", "AttendLL", timing$question)
  timing$question <- gsub("Attend-SS", "AttendSS", timing$question)
  timing <- timing |>
    tidyr::separate(question, c("index", "measure"), sep = "-") |>
    tidyr::spread(measure, value)
  timing$q[timing$index %in% c("I16")] <- 1
  timing$q[timing$index %in% c("I8", "I24")] <- 2
  timing$q[timing$index %in% c("I4", "I12", "I20", "I28")] <- 3
  timing$q[timing$index %in% c("I2", "I6", "I10", "I14", "I18", "I22", "I26", "I30")] <- 4
  timing$q[timing$index %in% c("I1", "I3", "I5", "I7", "I9", "I11", "I13", "I15", "I17", "I19",
                               "I21", "I23", "I25", "I27", "I29", "I31")] <- 5
  timing$q[timing$index %in% c("AttendSS", "AttendLL")] <- 6
  timing[4:7] <- sapply(timing[4:7], as.numeric)
  return(dplyr::arrange(timing, ResponseId, q))

}

#' Converts answers from 5.5 trial probability discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns
#'
#' @return A dataframe with the ResponseId, index, and response (sc or lu).
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' ans_pd(five.fivetrial_pd)
ans_pd <- function(df) {
  ans <- df |>
    dplyr::select(ResponseId, paste0("I", 1:31), dplyr::starts_with("Attend")) |>
    dplyr::select(-dplyr::contains("Timing"), -dplyr::contains("_DO")) %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "index", values_to = "response") %>%
    dplyr::filter(complete.cases(.)) |>
    dplyr::mutate(response = ifelse(stringr::str_detect(response, "for sure"), "sc", "lu"))
  ans$index <- gsub("-", "", ans$index)
  return(ans)
}

#' Calculate scores, answers, and timing for 5.5 trial probability discounting from Qualtrics template
#'
#' @param df A dataframe containing all the columns from the template.
#'
#' @return A dataframe with h/ep50 values, answers, timing
#' @export
#'
#' @examples
#' calc_pd(five.fivetrial_pd)
calc_pd <- function(df) {
  return(dplyr::left_join(timing_pd(df), ans_pd(df), by = c("ResponseId", "index")) |>
           dplyr::left_join(dplyr::select(score_pd(df),
                                          ResponseId, attentionflag,
                                          hval, etheta50, ep50),
                            by = c("ResponseId")) |>
           dplyr::arrange(ResponseId, q))
}
