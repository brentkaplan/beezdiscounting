## data-raw/mcq21_lookup.R
## Adds the internal `lookup21` table: the canonical Kirby & Marakovic (1996)
## Table 1 design for the 21-item MCQ, in the SAME ladder order used by the
## Kaplan et al. (2014) 21-item Automated Scorer ("Ordered Data" sheet).
## Row order is load-bearing: score_mcq()/inn()/prop_ss() walk this order.
## Note the ranks are not a clean 7x3 grid: rank 3 has no Large item and
## rank 4 has two — this matches both the paper and the Excel scorer.
## Run from the package root:  Rscript data-raw/mcq21_lookup.R

e <- new.env()
load("R/sysdata.rda", envir = e)
## `lookup` must exist; `lookup21` may exist from a prior run (rerunnable).
stopifnot("lookup" %in% ls(e), all(ls(e) %in% c("lookup", "lookup21")))
lookup <- get("lookup", envir = e)
stopifnot(nrow(lookup) == 27L)

lookup21 <- data.frame(
  questionid = c(
    4,
    15,
    7,
    20,
    9,
    12,
    8,
    16,
    14,
    10,
    3,
    18,
    11,
    2,
    19,
    21,
    6,
    17,
    5,
    13,
    1
  ),
  magnitude = c(
    "S",
    "M",
    "L",
    "S",
    "M",
    "L",
    "S",
    "M",
    "S",
    "M",
    "L",
    "L",
    "S",
    "M",
    "L",
    "S",
    "M",
    "L",
    "S",
    "M",
    "L"
  ),
  kindiff = c(
    0.0007,
    0.0007,
    0.0007,
    0.0032,
    0.0032,
    0.0031,
    0.0057,
    0.0055,
    0.0083,
    0.0089,
    0.0077,
    0.0086,
    0.0160,
    0.0150,
    0.0159,
    0.0250,
    0.0359,
    0.0375,
    0.1333,
    0.1292,
    0.1310
  ),
  k_rank = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7),
  ss_amount = c(
    34,
    53,
    83,
    27,
    48,
    65,
    21,
    47,
    30,
    40,
    67,
    50,
    25,
    40,
    45,
    16,
    32,
    40,
    15,
    24,
    30
  ),
  ll_amount = c(
    35,
    55,
    85,
    30,
    55,
    75,
    30,
    60,
    35,
    65,
    85,
    80,
    35,
    55,
    70,
    30,
    55,
    70,
    35,
    55,
    85
  ),
  delay = c(
    43,
    55,
    35,
    35,
    45,
    50,
    75,
    50,
    20,
    70,
    35,
    70,
    25,
    25,
    35,
    35,
    20,
    20,
    10,
    10,
    14
  ),
  stringsAsFactors = FALSE
)

## Integrity gates.
stopifnot(
  nrow(lookup21) == 21L,
  setequal(lookup21$questionid, 1:21),
  all(table(lookup21$magnitude) == 7L),
  all(lookup21$ll_amount > lookup21$ss_amount),
  ## per-magnitude ladder strictly ascending in table order
  all(vapply(
    split(seq_len(21), lookup21$magnitude),
    function(i) all(diff(lookup21$kindiff[i]) > 0),
    TRUE
  )),
  ## Kirby Eq. 1 within published rounding (2 sig figs -> loose gate)
  max(
    abs(
      (lookup21$ll_amount / lookup21$ss_amount - 1) /
        lookup21$delay -
        lookup21$kindiff
    ) /
      lookup21$kindiff
  ) <
    0.05
)

usethis::use_data(lookup, lookup21, internal = TRUE, overwrite = TRUE)

## Post-write verification.
e2 <- new.env()
load("R/sysdata.rda", envir = e2)
stopifnot(
  setequal(ls(e2), c("lookup", "lookup21")),
  identical(get("lookup", envir = e2), lookup),
  identical(get("lookup21", envir = e2), lookup21)
)
message("sysdata.rda rebuilt with lookup (27) + lookup21 (21)")
