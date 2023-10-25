## prepare 5.5 trial example data

five.fivetrial <- readr::read_csv("five.fivetrial.csv")
# remove 2 first two rows with header information
# remove SC0/SC1 scoring columns
five.fivetrial <- five.fivetrial[-c(1:2), -c(173, 174)]

usethis::use_data(five.fivetrial, overwrite = TRUE)
