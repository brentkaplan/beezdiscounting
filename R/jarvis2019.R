#' Jarvis (2019) delay-discounting indifference points
#'
#' Long-format delay-discounting indifference-point data from Jarvis et al.
#' (2019), used to demonstrate the scale-location-truncated beta (SLT-beta)
#' family's handling of boundary indifference points (exactly 0 or 1). Many
#' observations sit at the bounds, where ordinary beta regression is undefined.
#'
#' The original source has a non-unique participant identifier (78 unique values
#' across 126 records); here each record is treated as one subject-record and is
#' assigned a unique synthetic `id` (the row index of the source file). Use this
#' `id` as the subject grouping for [fit_dd_tmb()].
#'
#' @format A [tibble][tibble::tibble] with 882 rows (126 subjects x 7 delays)
#'   and 3 columns:
#' \describe{
#'   \item{id}{Factor; unique subject-record identifier.}
#'   \item{x}{Numeric; delay in days (7, 30, 180, 365, 730, 1460, 2920 =
#'     1 week through 8 years).}
#'   \item{y}{Numeric; indifference point as a proportion of the larger-later
#'     reward, in the closed interval `[0, 1]` (0 and 1 occur and are valid).}
#' }
#'
#' @source Jarvis, B. P., et al. (2019). Delay-discounting dataset
#'   (N = 126; 1 week - 8 years). Source CSV not redistributed with the package;
#'   see `data-raw/jarvis2019.R`.
#' @keywords datasets
"jarvis2019"
