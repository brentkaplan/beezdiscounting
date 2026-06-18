#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @importFrom generics augment
#' @export
generics::augment

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme ranef
#' @export
nlme::ranef

#' @importFrom beezdemand plot_qq
#' @export
beezdemand::plot_qq

utils::globalVariables(c(
    ".",
    "ResponseId",
    "response",
    "attentionflag",
    "kval",
    "ed50",
    "Attend-SS",
    "Attend-LL",
    "magnitude",
    "question",
    "measure",
    "value",
    "subjectid",
    "questionid",
    "hval",
    "etheta50",
    "ep50",
    "lookup",
    "k_rank",
    "Statistic",
    "overall_k",
    "composite_consistency",
    "geomean_k",
    "group",
    "metric",
    "estimate",
    "std.error",
    "conf",
    "model",
    "term",
    "fit",
    "tidy_summary",
    "glance_summary",
    "R2",
    "conf_int",
    "x",
    "y",
    "id",
    "method"
    ))
