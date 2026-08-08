## data-raw/mcq21_data.R
## Two illustrative 21-item MCQ subjects:
##   1: Kirby & Marakovic's (1996) median responder — SIR on ladder
##      positions 1-8, LL after (overall k = geomean(0.0055, 0.0083) ~ 0.0068)
##   2: uniformly shallow responder (all LL, overall k = 0.0007)
sir_qids <- c(4, 15, 7, 20, 9, 12, 8, 16)
mcq21 <- data.frame(
  subjectid = rep(1:2, each = 21),
  questionid = rep(1:21, times = 2),
  response = c(as.integer(!(1:21 %in% sir_qids)), rep(1L, 21))
)
usethis::use_data(mcq21, overwrite = TRUE)
