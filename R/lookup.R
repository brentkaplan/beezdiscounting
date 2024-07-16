#' MCQ Lookup Table
#'
#' This dataset contains a lookup table for the 27-item Monetary Choice
#' Questionnaire (MCQ), which assesses delay discounting behavior.
#' Each row represents a question with associated metadata.
#'
#' @format A data frame with 27 rows and 4 variables:
#' \describe{
#'   \item{questionid}{Integer. Unique identifier for each question.}
#'   \item{magnitude}{Character. Indicates the magnitude of the reward, either "S" (small), "M" (medium), or "L" (large).}
#'   \item{kindiff}{Numeric. The difference in subjective value between the immediate and delayed rewards.}
#'   \item{k_rank}{Numeric. The rank-based k-value, representing the degree of delay discounting.}
#' }
#' @source See Kaplan et al. (2016) for more information.
"lookup"
