## Extract the item table of the Kaplan et al. (2014) 21-item MCQ Excel scorer
## as a test fixture (audit F-BZ1-1). Source (not bundled):
##   ~/Dropbox/sharedDerek/Automating_Kirby_Discounting/Scorer-CurrentVersions/
##   21-Item-MCQ-Auto-Scorer-Kaplan_et_al-2014-20180419.xlsx
## Sheet "Ordered Data" A3:C23 = ladder order (trial #, magnitude, k at
## indifference, stored at 4 dp); sheet "Enter Data Here" A4:B24 = each
## trial's delay ("N days from now"). Sheet "All" B24 = 0.1333 is the overall
## ladder edge, and All!C54 assigns that edge itself (not a geometric mean) to
## the all-smaller-sooner switch point. The workbook stores no reward amounts.
f <- path.expand(file.path(
  "~/Dropbox/sharedDerek/Automating_Kirby_Discounting/Scorer-CurrentVersions",
  "21-Item-MCQ-Auto-Scorer-Kaplan_et_al-2014-20180419.xlsx"
))
ord <- readxl::read_excel(f, sheet = "Ordered Data", range = "A3:C23",
                          col_names = c("questionid", "magnitude", "kindiff"))
dly <- readxl::read_excel(f, sheet = "Enter Data Here", range = "A4:B24",
                          col_names = c("questionid", "text"))
dly$delay <- as.numeric(sub(" days from now", "", dly$text))
out <- data.frame(
  ladder_pos = seq_len(nrow(ord)),
  questionid = as.integer(ord$questionid),
  magnitude = ord$magnitude,
  kindiff = ord$kindiff,
  delay = dly$delay[match(ord$questionid, as.numeric(dly$questionid))]
)
edge <- as.numeric(readxl::read_excel(f, sheet = "All", range = "B24",
                                      col_names = FALSE)[[1]])
stopifnot(nrow(out) == 21L, !anyNA(out), identical(edge, 0.1333))
utils::write.csv(out, "tests/testthat/fixtures/mcq21/kaplan2014-21item-ordered.csv",
                 row.names = FALSE)
