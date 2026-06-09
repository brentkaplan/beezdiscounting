## data-raw/mcq27_lookup.R
## Completes the internal `lookup` table with the canonical Kirby, Petry, &
## Bickel (1999) Table 3 stimulus design (ss_amount / ll_amount / delay).
## The mapping is ORDER-PRESERVING: `lookup` ships in k-rank order
## (13, 1, 9, 20, ...), and score_mcq27()/inn()/prop_ss() depend on that order,
## so we index the design into lookup by questionid rather than merge() (base
## merge re-sorts by the join key and would corrupt scoring).
## Run from the package root:  Rscript data-raw/mcq27_lookup.R

## 1. Load the current internal data; assert `lookup` is the only object so the
##    use_data() re-save (which overwrites all of R/sysdata.rda) drops nothing.
e <- new.env()
load("R/sysdata.rda", envir = e)
stopifnot(setequal(ls(e), "lookup"))
lookup <- get("lookup", envir = e)
base_cols <- c("questionid", "magnitude", "kindiff", "k_rank")
stopifnot(all(base_cols %in% names(lookup)), nrow(lookup) == 27L)
## Rerunnable: reduce to the canonical base columns, dropping any previously
## added ss_amount/ll_amount/delay, so re-running on a completed lookup is a no-op.
lookup <- lookup[, base_cols]
old <- lookup

## 2. Canonical Kirby (1999) Table 3 design: questionid -> SIR / LDR / delay(days).
##    Cross-verified against Kaplan et al. (2016) MCQ Scorer (Fig 2 + Table 1).
design <- data.frame(
  questionid = c(13, 1, 9, 20, 6, 17, 26, 24, 12, 22, 16, 15, 3, 10, 2,
                 18, 21, 25, 5, 14, 23, 7, 8, 19, 11, 27, 4),
  ss_amount  = c(34, 54, 78, 28, 47, 80, 22, 54, 67, 25, 49, 69, 19, 40, 55,
                 24, 34, 54, 14, 27, 41, 15, 25, 33, 11, 20, 31),
  ll_amount  = c(35, 55, 80, 30, 50, 85, 25, 60, 75, 30, 60, 85, 25, 55, 75,
                 35, 50, 80, 25, 50, 75, 35, 60, 80, 30, 55, 85),
  delay      = c(186, 117, 162, 179, 160, 157, 136, 111, 119, 80, 89, 91, 53, 62, 61,
                 29, 30, 30, 19, 21, 20, 13, 14, 14, 7, 7, 7),
  stringsAsFactors = FALSE
)

## 3. Integrity gates BEFORE touching lookup.
stopifnot(
  !anyDuplicated(design$questionid),          # no duplicate items
  setequal(design$questionid, old$questionid) # exact 27-id coverage, no orphans
)

## 4. Map the design into the EXISTING lookup row order (do NOT merge()).
i <- match(lookup$questionid, design$questionid)
lookup$ss_amount <- design$ss_amount[i]
lookup$ll_amount <- design$ll_amount[i]
lookup$delay     <- design$delay[i]

## 5. Post-map gates.
stopifnot(
  !anyNA(lookup$ss_amount), !anyNA(lookup$ll_amount), !anyNA(lookup$delay),
  identical(lookup$questionid, old$questionid), # row order preserved
  all(lookup$ll_amount > lookup$ss_amount)      # LDR always exceeds SIR
)

## Kirby Eq. 1 transcription gate: k = (LDR/SIR - 1) / delay must reproduce kindiff.
k_eq1   <- (lookup$ll_amount / lookup$ss_amount - 1) / lookup$delay
rel_dif <- abs(k_eq1 - lookup$kindiff) / lookup$kindiff
if (max(rel_dif) > 1e-5) {
  bad <- lookup$questionid[rel_dif > 1e-5]
  stop("mcq27 design inconsistent with lookup$kindiff (Kirby Eq.1) for questionid: ",
       paste(bad, collapse = ", "), call. = FALSE)
}
message(sprintf("Eq.1 check OK: max relative diff = %.2e", max(rel_dif)))

## 6. Re-save sysdata.rda (use_data overwrites the whole file; lookup is its only object).
usethis::use_data(lookup, internal = TRUE, overwrite = TRUE)

## 7. Post-write verification: reload from disk and re-assert.
e2 <- new.env()
load("R/sysdata.rda", envir = e2)
chk <- get("lookup", envir = e2)
stopifnot(
  setequal(ls(e2), "lookup"),
  ncol(chk) == 7L,
  identical(chk$questionid, old$questionid),
  all(c("ss_amount", "ll_amount", "delay") %in% names(chk))
)
message("sysdata.rda rebuilt: lookup is now ", nrow(chk), " x ", ncol(chk),
        " (", paste(names(chk), collapse = ", "), ")")
