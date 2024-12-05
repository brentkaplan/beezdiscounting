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

## make mock discounting data
set.seed(123)
num_participants <- 100
delay_values <- c(1, 7, 30, 90, 180, 365)

dd_ip <- do.call(rbind, lapply(1:num_participants, function(participant) {
  k <- runif(1, 0.001, 1) # random discount rate for each participant
  y <- 1 / (1 + k * delay_values) + rnorm(length(delay_values), 0, 0.05) # Indifference points
  data.frame(
    id = paste0("P", participant),
    x = delay_values,
    y = pmax(0, y) # indifference points are non-negative
  )
}))

usethis::use_data(dd_ip, overwrite = TRUE)
