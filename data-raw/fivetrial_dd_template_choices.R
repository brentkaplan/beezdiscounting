## Extract every choice of the bundled 5.5-trial delay-discounting template
## (inst/5.5_Trial_Discounting_Template_1k.qsf) with its exported numeric code
## (RecodeValues) and its exported text (HTML tags stripped, as Qualtrics does).
## Ground truth for the normalize_dd_response() tests (audit F-BZ2-3 / F-BZ2-2).
qsf <- jsonlite::fromJSON("inst/5.5_Trial_Discounting_Template_1k.qsf",
                          simplifyVector = FALSE)
rows <- list()
for (el in qsf$SurveyElements) {
  if (!identical(el$Element, "SQ")) next
  p <- el$Payload
  tag <- p$DataExportTag
  if (!grepl("^(I[0-9]+|Attend-(SS|LL))$", tag)) next
  for (cid in names(p$Choices)) {
    rows[[length(rows) + 1L]] <- data.frame(
      item = tag,
      code = p$RecodeValues[[cid]] %||% cid,
      text = gsub("<[^>]+>", "", p$Choices[[cid]]$Display),
      # SS = the immediate option ("... now") in every item, incl. both
      # attention checks (Attend-SS: $0 now; Attend-LL: the full amount now)
      choice = if (grepl("now", p$Choices[[cid]]$Display)) "ss" else "ll"
    )
  }
}
out <- do.call(rbind, rows)
utils::write.csv(out, "tests/testthat/fixtures/fivetrial/dd_template_choices.csv",
                 row.names = FALSE)
