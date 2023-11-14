## prepare 5.5 trial dd example data

five.fivetrial_dd <- readr::read_csv("five.fivetrial_dd.csv")
# remove 2 first two rows with header information
# remove SC0/SC1 scoring columns
five.fivetrial_dd <- five.fivetrial_dd[-c(1:2), -c(173, 174)]
five.fivetrial_dd$ResponseId <- 1:4

usethis::use_data(five.fivetrial_dd, overwrite = TRUE)

## prepare 5.5 trial pd example data
five.fivetrial_pd <- readr::read_csv("five.fivetrial_pd.csv")
five.fivetrial_pd <- five.fivetrial_pd[five.fivetrial_pd$ResponseId %in%
                                         c(sample(five.fivetrial_pd$ResponseId, 3),
                                           "R_20MSveDnlaykgYW"), ]
five.fivetrial_pd$ResponseId <- 1:4

usethis::use_data(five.fivetrial_pd, overwrite = TRUE)
