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

#' Check for Unsystematic Data Violations
#'
#' This function checks a dataset for violations of two criteria commonly used to identify unsystematic delay-discounting data:
#' - Criterion 1: Any subsequent value of `y` exceeds the previous value by more than a specified proportion of the larger later reward (`ll`).
#' - Criterion 2: The last value of `y` is not at least a specified proportion less than the first value of `y`.
#'
#' @param dat A data frame containing the delay-discounting data. It must have at least two columns:
#'   - `id`: A unique identifier for the data set.
#'   - `y`: The indifference points to be analyzed.
#' @param ll A numeric value representing the larger later reward. Default is 1.
#' @param c1 A numeric value for the threshold proportion for Criterion 1. Default is 0.2.
#' @param c2 A numeric value for the threshold proportion for Criterion 2. Default is 0.1.
#'
#' @return A tibble with the following columns:
#'   - `id`: The unique identifier for the data set.
#'   - `c1_violation`: Logical value indicating whether Criterion 1 was violated.
#'   - `c2_violation`: Logical value indicating whether Criterion 2 was violated.
#' @export
#'
#' @examples
#' data <- tibble::tibble(
#'   id = c(rep("P1", 6)),
#'   x = c(1, 7, 30, 90, 180, 365), # delays
#'   y = c(0.9, 0.5, 0.3, 0.2, 0.1, 0.05) # indifference points
#' )
#' check_unsystematic(data, ll = 1, c1 = 0.2, c2 = 0.1)
check_unsystematic <- function(dat, ll = 1, c1 = .2, c2 = .1) {

  c1_threshold <- c1 * ll
  c2_threshold <- c2 * ll

  # C1: Check if any subsequent y exceeds the previous y by more than c1 of ll
  c1_check <- any(diff(dat$y) > c1_threshold / ll)

  # C2: Check if the last y is not at least c2 less than the first y
  c2_check <- (dat$y[nrow(dat)] >= (dat$y[1] - c2_threshold / ll))

  out <- tibble::tibble(
    id = unique(dat$id),
    c1_violation = c1_check,
    c2_violation = c2_check
  )

  return(out)

}


#' Calculate R-Squared for a Model
#'
#' This function calculates the coefficient of determination (\eqn{R^2}) for a given model by comparing the sum of squared errors (SSE)
#' to the total sum of squares (SST).
#'
#' @param model A fitted model object. The model must have `resid()` and `fitted()` methods to extract residuals and fitted values.
#'
#' @return A numeric value representing the \eqn{R^2} value of the model. Returns `NA` if the model is `NULL`.
#' @importFrom stats resid fitted
#' @export
#'
#' @examples
#' # Example using a simple linear model
#' data <- data.frame(x = 1:10, y = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11))
#' lm_model <- lm(y ~ x, data = data)
#' calc_r2(lm_model)
calc_r2 <- function(model) {
  if (!is.null(model)) {
    residuals <- resid(model)
    fitted_values <- fitted(model)
    observed_values <- fitted_values + residuals
    sse <- sum(residuals^2)
    sst <- sum((observed_values - mean(observed_values))^2) # Total sum of squares
    r_squared <- 1 - (sse / sst)
    return(r_squared)
  } else {
    return(NA) # Handle cases where the model is NULL
  }
}


#' Calculate Confidence Intervals for a Parameter
#'
#' This function computes the lower and upper bounds of the confidence interval for a parameter estimate, given its standard error,
#' a specified significance level, and the degrees of freedom from the model.
#'
#' @param estimate A numeric value representing the parameter estimate.
#' @param std_error A numeric value representing the standard error of the parameter estimate.
#' @param model A fitted model object that provides the residual degrees of freedom via `df.residual()`.
#' @param alpha A numeric value representing the significance level. Default is 0.05 (95% confidence interval).
#'
#' @return A numeric vector of length two:
#'   - First element: Lower bound of the confidence interval.
#'   - Second element: Upper bound of the confidence interval.
#' @importFrom stats qt df.residual
#' @export
#'
#' @examples
#' # Example using a linear model
#' data <- data.frame(x = 1:10, y = c(2.3, 2.1, 3.7, 4.5, 5.1, 6.8, 7.3, 7.9, 9.2, 10.1))
#' lm_model <- lm(y ~ x, data = data)
#' calc_conf_int(estimate = 0.5, std_error = 0.1, model = lm_model, alpha = 0.05)
calc_conf_int <- function(estimate, std_error, model, alpha = 0.05) {
  t_crit <- qt(1 - alpha / 2, df = df.residual(model)) # critical t-value
  lower <- estimate - t_crit * std_error
  upper <- estimate + t_crit * std_error
  c(lower, upper)
}


#' Calculate Area-Under-the-Curve (AUC) Metrics for Delay Discounting Data
#'
#' This function calculates three types of Area-Under-the-Curve (AUC) metrics for delay discounting data:
#' regular AUC (using raw delays), log10 AUC (using logarithmically scaled delays), and rank AUC (using ordinally scaled delays).
#' These metrics provide different perspectives on the rate of delay discounting.
#'
#' @param dat A data frame containing delay discounting data.
#'   It must include the following columns:
#'   - `id`: Participant or group identifier.
#'   - `x`: Delay values (e.g., in days).
#'   - `y`: Indifference point values (e.g., subjective value of the delayed reward).
#'
#' @return A tibble with the following columns:
#'   - `id`: The participant or group identifier.
#'   - `auc_regular`: The regular AUC, calculated using the raw delay values.
#'   - `auc_log10`: The log10 AUC, calculated using logarithmically transformed delay values.
#'   - `auc_rank`: The rank AUC, calculated using ordinally scaled delay values.
#'
#' @export
#'
#' @examples
#' # Example data
#' data <- data.frame(
#'   id = rep("P1", 6),
#'   x = c(1, 7, 30, 90, 180, 365),
#'   y = c(0.8, 0.5, 0.3, 0.2, 0.1, 0.05)
#' )
#'
#' # Calculate AUC metrics for a single participant
#' calc_aucs(data)
calc_aucs <- function(dat) {
  dat <- dat |> dplyr::arrange(x)  # Ensure dat is sorted by delay (x)

  # Regular AUC
  auc_regular_raw <- sum(
    (diff(dat$x) * (dat$y[-length(dat$y)] + dat$y[-1]) / 2)
  )
  max_possible_area <- diff(range(dat$x))  # Total range of x
  auc_regular <- auc_regular_raw / max_possible_area  # Normalize

  # Log10 AUC
  log_x <- log10(dat$x + 1)  # Add 1 to handle log10(0) issue
  log_x <- log_x / max(log_x)  # Scale to proportions
  auc_log <- sum(
    (diff(log_x) * (dat$y[-length(dat$y)] + dat$y[-1]) / 2)
  )

  # Rank AUC
  rank_x <- seq_len(nrow(dat)) / nrow(dat)  # Ranks as proportions
  auc_rank <- sum(
    (diff(rank_x) * (dat$y[-length(dat$y)] + dat$y[-1]) / 2)
  )

  # Return a tibble with the results
  tibble::tibble(
    id = unique(dat$id),
    auc_regular = auc_regular,
    auc_log10 = auc_log,
    auc_rank = auc_rank
  )
}
