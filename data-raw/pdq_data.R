## data-raw/pdq_data.R
## Two illustrative 30-item PDQ subjects:
##   1: consistent moderate discounter -- guaranteed through h rank 5,
##      risky after, in every block (block h = sqrt(h5 * h6) ~ 1.22)
##   2: uniformly risk-tolerant responder (all risky; block h = ladder
##      minima, the smallest h consistent with the pattern)
pdq <- data.frame(
  subjectid = rep(1:2, each = 30),
  questionid = rep(1:30, times = 2),
  response = c(rep(rep(c(0L, 1L), each = 5), times = 3), rep(1L, 30))
)
usethis::use_data(pdq, overwrite = TRUE)
