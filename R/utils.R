#' Generate fake MCQ data
#'
#' @param n_ids Number of subjectids
#' @param n_items Number of trials
#' @param seed Random seed
#' @param prop.na Proportion of NAs in the entire data set
#'
#' @return Dataframe of subjectid, questionid, and response
#' @export
#'
#' @examples
#' generate_data_mcq(n_ids = 2, n_items = 27, prop.na = .01)
generate_data_mcq <- function(n_ids = 100, n_items = 27, seed = 1234, prop.na = 0) {
  set.seed(seed)
  fake_data <- data.frame(
    subjectid = rep(1:n_ids, each = n_items),
    questionid = rep(1:n_items, times = n_ids),
    response = sample(0:1, n_ids*n_items, replace = TRUE)
  )

  fake_data$response[sample(nrow(fake_data), round(nrow(fake_data)*prop.na))] <- NA
  return(fake_data)

}
