#' Parameter naming registry for beezdiscounting mixed-effects models
#'
#' Single source of truth for the canonical parameters of the IP-family
#' mixed-effects discounting model: the subject discount rate `k`, the
#' SLT-beta precision `phi`, and the hyperboloid nonlinearity exponent `s`
#' (Green-Myerson / Rachlin). Fields mirror `beezdemand`'s
#' `.beezdemand_param_registry` (canonical, description, constraint,
#' valid_scales, default_scale).
#'
#' @keywords internal
#' @noRd
.dd_param_registry <- list(
  k = list(
    canonical = "k",
    description = "Discount rate (Mazur hyperbolic / exponential)",
    constraint = "k > 0",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  ),

  phi = list(
    canonical = "phi",
    description = "SLT-beta precision (population)",
    constraint = "phi > 0",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  ),

  s = list(
    canonical = "s",
    description = "Hyperboloid nonlinearity exponent (Green-Myerson / Rachlin)",
    constraint = "s > 0",
    valid_scales = c("natural", "log", "log10"),
    default_scale = "natural"
  )
)
