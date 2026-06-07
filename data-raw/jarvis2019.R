## ----------------------------------------------------------------------------
## Build the bundled, de-duplicated Jarvis (2019) delay-discounting fixture.
##
## SOURCE: Data_brent.csv (Jarvis et al. 2019; N = 126, delays 1wk..8yr).
## The CSV is NOT shipped in the package (private source data); it lives at
## ~/Downloads/Data_brent.csv. Run this script manually to (re)build
## data/jarvis2019.rda. PID is NON-UNIQUE (78 unique / 126 rows), so we assign
## a synthetic unique subject id (the row index) and treat each row as one
## subject-record. A mixed model that grouped duplicated PIDs as one subject
## would mis-group; the unique `id` here prevents that.
## ----------------------------------------------------------------------------

csv_path <- path.expand("~/Downloads/Data_brent.csv")
if (!file.exists(csv_path)) {
  stop(
    "Source CSV not found at ", csv_path, ".\n",
    "Place the Jarvis (2019) Data_brent.csv there and re-run this script.\n",
    "The CSV is not distributed with the package.",
    call. = FALSE
  )
}

d <- utils::read.csv(csv_path, check.names = FALSE)
delays <- c(7, 30, 180, 365, 730, 1460, 2920)  # 1wk,1mo,6mo,1yr,2yr,4yr,8yr (days)

# Each row is one subject-record; PID is non-unique so we use the row index.
n <- nrow(d)
ip_cols <- d[, -1, drop = FALSE]               # drop the PID column
stopifnot(ncol(ip_cols) == length(delays))

jarvis2019 <- tibble::tibble(
  id = factor(rep(seq_len(n), each = length(delays))),
  x  = rep(delays, times = n),
  y  = as.numeric(t(as.matrix(ip_cols)))
)

# Sanity: boundary IPs are the point of this fixture.
message(sprintf(
  "jarvis2019: %d subjects, %d obs, IP in [%.3f, %.3f]; #(y==0)=%d, #(y==1)=%d",
  n, nrow(jarvis2019), min(jarvis2019$y), max(jarvis2019$y),
  sum(jarvis2019$y == 0), sum(jarvis2019$y == 1)
))

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(jarvis2019, overwrite = TRUE)
} else {
  # Fallback: save directly; run this script from the package root so that
  # the relative path "data/jarvis2019.rda" resolves correctly.
  if (!dir.exists("data")) dir.create("data", recursive = TRUE)
  save(jarvis2019, file = "data/jarvis2019.rda")
  message("Saved to data/jarvis2019.rda")
}
