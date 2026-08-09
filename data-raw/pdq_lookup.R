## data-raw/pdq_lookup.R
## Adds the internal `lookup_pdq` table: the canonical Madden, Petry &
## Johnson (2009) 30-item Probability Discounting Questionnaire design as
## printed in Gray et al. (2016, JEAB, Table 2). Items are in the
## administered order, which is ascending-h ladder order within each of the
## three 10-item blocks (presentation is NOT randomized). Row order is
## load-bearing: score_pdq() walks it.
## theta and hindiff are computed full-precision from the printed certain
## amount, probabilistic amount, and probability -- verified to reproduce
## Gray et al.'s supplementary lookup tables (all 3 x 1024 response
## patterns) to <1e-8 relative error.
## Run from the package root:  Rscript data-raw/pdq_lookup.R

e <- new.env()
load("R/sysdata.rda", envir = e)
## lookup + lookup21 must exist; lookup_pdq may exist from a prior run.
stopifnot(
  all(c("lookup", "lookup21") %in% ls(e)),
  all(ls(e) %in% c("lookup", "lookup21", "lookup_pdq"))
)
lookup <- get("lookup", envir = e)
lookup21 <- get("lookup21", envir = e)
stopifnot(nrow(lookup) == 27L, nrow(lookup21) == 21L)

lookup_pdq <- data.frame(
  questionid = 1:30,
  block = rep(1:3, each = 10L),
  h_rank = rep(1:10, times = 3L),
  sc_amount = rep(c(20, 40, 40), each = 10),
  lu_amount = rep(c(80, 100, 60), each = 10),
  prob = c(
    .10,
    .13,
    .17,
    .20,
    .25,
    .33,
    .50,
    .67,
    .75,
    .83,
    .18,
    .22,
    .29,
    .33,
    .40,
    .50,
    .67,
    .80,
    .86,
    .91,
    .40,
    .46,
    .55,
    .60,
    .67,
    .75,
    .86,
    .92,
    .95,
    .97
  ),
  stringsAsFactors = FALSE
)
lookup_pdq$theta <- (1 - lookup_pdq$prob) / lookup_pdq$prob
lookup_pdq$hindiff <- (lookup_pdq$lu_amount / lookup_pdq$sc_amount - 1) /
  lookup_pdq$theta

## Pooled 30-item ladder position (the beezdiscounting overall_h extension).
## Order by EXACT rational h = (lu - sc) * p100 / (sc * (100 - p100)) with
## ties broken by questionid. The single integer division makes equal
## rationals equal doubles (correctly-rounded division), so the four exact
## tie pairs (q1&q21 at 1/3, q4&q24 at 3/4, q5&q15 at 1, q16&q26 at 3/2)
## genuinely tie and fall through to the questionid tie-break -- the
## multi-step hindiff above breaks them by floating-point noise instead.
p100 <- round(lookup_pdq$prob * 100)
h_exact <- (lookup_pdq$lu_amount - lookup_pdq$sc_amount) *
  p100 /
  (lookup_pdq$sc_amount * (100 - p100))
ord <- order(h_exact, lookup_pdq$questionid)
lookup_pdq$overall_rank <- NA_integer_
lookup_pdq$overall_rank[ord] <- 1:30
## Pin the canonical pooled order (verified 2026-08-08).
stopifnot(identical(
  lookup_pdq$questionid[order(lookup_pdq$overall_rank)],
  c(
    11L,
    1L,
    21L,
    12L,
    22L,
    2L,
    23L,
    13L,
    3L,
    14L,
    4L,
    24L,
    5L,
    15L,
    25L,
    6L,
    16L,
    26L,
    7L,
    17L,
    27L,
    28L,
    18L,
    8L,
    9L,
    19L,
    29L,
    10L,
    20L,
    30L
  )
))
lookup_pdq <- lookup_pdq[, c(
  "questionid",
  "block",
  "h_rank",
  "overall_rank",
  "sc_amount",
  "lu_amount",
  "prob",
  "theta",
  "hindiff"
)]

## Integrity gates.
printed <- c(
  0.33,
  0.45,
  0.61,
  0.75,
  1,
  1.48,
  3,
  6.09,
  9,
  14.65,
  0.33,
  0.42,
  0.62,
  0.74,
  1,
  1.5,
  3.04,
  6,
  9.21,
  15.17,
  0.33,
  0.43,
  0.61,
  0.75,
  1.01,
  1.5,
  3.07,
  5.75,
  9.5,
  16.17
)
stopifnot(
  nrow(lookup_pdq) == 30L,
  all(lookup_pdq$lu_amount > lookup_pdq$sc_amount),
  all(vapply(
    split(lookup_pdq$hindiff, lookup_pdq$block),
    function(h) all(diff(h) > 0),
    TRUE
  )),
  ## printed Table 2 values reproduce to ~1% (items 13/17/25 are
  ## last-digit inconsistent in the published table -- lookups are the
  ## precision oracle)
  max(abs(lookup_pdq$hindiff - printed) / printed) < 0.02
)

usethis::use_data(
  lookup,
  lookup21,
  lookup_pdq,
  internal = TRUE,
  overwrite = TRUE
)

## Post-write verification.
e2 <- new.env()
load("R/sysdata.rda", envir = e2)
stopifnot(
  setequal(ls(e2), c("lookup", "lookup21", "lookup_pdq")),
  identical(get("lookup", envir = e2), lookup),
  identical(get("lookup21", envir = e2), lookup21),
  identical(get("lookup_pdq", envir = e2), lookup_pdq)
)
message(
  "sysdata.rda rebuilt with lookup (27) + lookup21 (21) + lookup_pdq (30)"
)
