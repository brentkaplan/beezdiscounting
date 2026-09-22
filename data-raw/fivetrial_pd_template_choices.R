## Extract every choice of the 5.5-trial probability-discounting Qualtrics
## template with its exported numeric code (RecodeValues) and its exported text
## (HTML tags stripped, as Qualtrics does). Ground truth for the
## normalize_pd_response() tests (audit F-BZ2-3 / F-BZ2-2, PD follow-up).
##
## Source (not bundled): 55_Trial_Discounting_Probability_Template_100.qsf.
## Two copies were checked on 2026-09-22 and give identical choice tables:
##   ~/Dropbox/Work/VT/VTCRI/Qualtrics/reference-surveys/
##   ~/Dropbox/Work/UK/ (..._Template_100_-_Share.qsf)
## Every I-item and Attend-SS codes 1 = "... for sure" (sc), 2 = "... chance"
## (lu); Attend-LL lists the 1%-chance option first (1 = lu, 2 = sc).
`%||%` <- function(a, b) if (is.null(a)) b else a
qsf_path <- path.expand(file.path(
  "~/Dropbox/Work/VT/VTCRI/Qualtrics/reference-surveys",
  "55_Trial_Discounting_Probability_Template_100.qsf"
))
qsf <- jsonlite::fromJSON(qsf_path, simplifyVector = FALSE)
rows <- list()
for (el in qsf$SurveyElements) {
  if (!identical(el$Element, "SQ")) next
  p <- el$Payload
  tag <- p$DataExportTag
  if (!grepl("^(I[0-9]+|Attend-(SS|LL))$", tag)) next
  for (cid in names(p$Choices)) {
    txt <- gsub("<[^>]+>", "", p$Choices[[cid]]$Display)
    rows[[length(rows) + 1L]] <- data.frame(
      item = tag,
      code = p$RecodeValues[[cid]] %||% cid,
      text = txt,
      # SC = the certain option ("... for sure") in every item, incl. both
      # attention checks (Attend-SS: $0 for sure; Attend-LL: the full amount
      # for sure vs. a 1% chance)
      choice = if (grepl("for sure", txt)) "sc" else "lu"
    )
  }
}
out <- do.call(rbind, rows)
stopifnot(nrow(out) == 66L, all(table(out$item) == 2L))
utils::write.csv(out, "tests/testthat/fixtures/fivetrial/pd_template_choices.csv",
                 row.names = FALSE)
