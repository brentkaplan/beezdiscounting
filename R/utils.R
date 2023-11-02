#' Generate fake MCQ data
#'
#' @param n_ids Number of subjectids
#' @param n_items Number of trials
#' @param seed Random seed
#' @param prop_na Proportion of NAs in the entire data set
#'
#' @return Dataframe of subjectid, questionid, and response
#' @export
#'
#' @examples
#' generate_data_mcq(n_ids = 2, n_items = 27, prop_na = .01)
generate_data_mcq <- function(n_ids = 100, n_items = 27,
                              seed = 1234, prop_na = 0) {
  set.seed(seed)
  fake_data <- data.frame(
    subjectid = rep(1:n_ids, each = n_items),
    questionid = rep(1:n_items, times = n_ids),
    response = sample(0:1, n_ids * n_items, replace = TRUE)
  )

  fake_data$response[sample(nrow(fake_data),
                            round(nrow(fake_data) * prop_na))] <- NA
  return(fake_data)

}

#' Reshape MCQ data wide to long
#'
#' @param dat Wide format MCQ assuming subject id is in column 1
#' @param items Number of MCQ questions
#'
#' @return Long format data frame
#' @export
#'
wide_to_long_mcq <- function(dat, items = 27) {
  return(tidyr::pivot_longer(dat, cols = 2:ncol(dat),
                             names_to = "questionid",
                             values_to = "response") |>
    dplyr::mutate(questionid = rep(1:items,
                                   times = length(unique(dat$subjectid)))))
}



#' Reshape MCQ data long to wide
#'
#' @param dat Long format MCQ
#' @param q_col Name of the question column (default is "questionid")
#' @param ans_col Name of the answer column (defualt is "response")
#'
#' @return Wide format data frame
#' @export
#'
long_to_wide_mcq <- function(dat, q_col = "questionid", ans_col = "response") {
  return(tidyr::pivot_wider(dat, names_from = tidyr::all_of(q_col),
                            values_from = tidyr::all_of(ans_col)))
}


#' Reshape MCQ data from wide (as used in the 21- and 27-Item
#' Monetary Choice Questionnaire Automated Scorer) to long
#'
#' @param dat Wide format MCQ data as used in the Excel Automated Scorers
#'
#' @return Long format data frame
#' @export
#'
#' @examples
#' wide_to_long_mcq_excel(long_to_wide_mcq_excel(generate_data_mcq(2)))
#'
wide_to_long_mcq_excel <- function(dat) {
  return(tidyr::pivot_longer(dat, cols = 2:ncol(dat),
                      names_to = "subjectid",
                      values_to = "response") |>
         dplyr::select(subjectid, questionid, response) |>
         dplyr::arrange(subjectid))
}

#' Reshape MCQ data from long to wide (as used in the
#' 21- and 27-Item Monetary Choice Questionnaire Automated Scorer)
#'
#' @param dat Long format MCQ data
#' @param subj_col Character column name of subject ids
#' @param ans_col Character column name of responses
#'
#' @return Wide format MCQ data that can be used in the Excel Automated Scorers
#' @export
#'
#' @examples
#' long_to_wide_mcq_excel(generate_data_mcq(2))
long_to_wide_mcq_excel <- function(dat, subj_col = "subjectid",
                                   ans_col = "response") {
    return(tidyr::pivot_wider(dat, names_from = tidyr::all_of(subj_col),
                              values_from = tidyr::all_of(ans_col)))
}
