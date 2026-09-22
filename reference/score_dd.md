# Score 5.5 trial delay discounting from Qualtrics template

Score 5.5 trial delay discounting from Qualtrics template

## Usage

``` r
score_dd(df)
```

## Arguments

- df:

  A dataframe containing all the columns

## Value

A dataframe with id, indexes, response, k value, and effective delay 50.

## Details

Currently assumes the attending questions are present and labeled
"Attend-LL" and "Attend-SS"

Responses may be exported from Qualtrics either as choice text or as the
template's numeric codes. Numeric codes are read per item: every `I`
item and `Attend-SS` code 1 = the immediate ("now") option and 2 = the
delayed option, while `Attend-LL` lists the delayed option first (1 =
"in 25 years", 2 = "now"). Blank cells are treated as unanswered and
dropped; any other unrecognised value is an error.

`kval` is in units of 1/day and `ed50 = 1 / kval` in days. The terminal
k for each item is the reciprocal of the ED50 implied by the template's
31-delay tree, with a year of 365.25 days. The month-based delays (items
I17 to I23) were computed with a month of about 30.44 days rather than
365.25 / 12 = 30.4375, so those k values differ from the exact-month
values by a relative 8e-5 (well below the precision of any published
5.5-trial analysis). The template values are kept unchanged.

## Examples

``` r
score_dd(five.fivetrial_dd)
#> # A tibble: 4 × 6
#>   ResponseId index    response     kval attentionflag    ed50
#>        <int> <chr>    <chr>       <dbl> <chr>           <dbl>
#> 1          1 I21      ss        0.00671 No            149.   
#> 2          2 I5       ss        4.90    No              0.204
#> 3          3 AttendSS ss       NA       Yes            NA    
#> 4          4 AttendLL ll       NA       Yes            NA    
```
