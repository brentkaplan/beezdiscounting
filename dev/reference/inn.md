# Calculates item nearest neighbor imputation approach discussed by Yeh et al. (2023)

Calculates item nearest neighbor imputation approach discussed by Yeh et
al. (2023)

## Usage

``` r
inn(dat, random, verbose)
```

## Arguments

- dat:

  A single subject's 27-item MCQ data in long form

- random:

  Boolean whether to insert a random draw (0 or 1) for NAs

- verbose:

  Boolean whether to print subject and question ids pertaining to
  missing data

## Value

An imputed data set to be scored
