# Scoring the 5.5-trial discounting task

## A one-minute discounting task

The 5-trial adjusting-delay task (Koffarnus & Bickel, 2014) estimates a
discount rate in under a minute. Rather than titrating an amount at each
of many delays, it presents five binary choices in which the delay of
the larger reward adjusts based on the previous response, converging on
the delay at which the person is indifferent. The “5.5-trial” family
extends the original with a sixth attention-check trial and adds a
probability-discounting variant and a wider range of delays (Koffarnus,
Rzeszutek & Kaplan, 2021). All of it is delivered through the Qualtrics
minute discounting template (Koffarnus, Kaplan & Stein, 2017), and
`beezdiscounting` scores the template’s raw export directly.

## The raw data

The package includes one example export per task. They are wide
Qualtrics files: one row per respondent and many columns (the choice
items `I1`, `I3`, …, the per-item response-timing columns, and the
attention-check items).

``` r

data(five.fivetrial_dd)
dim(five.fivetrial_dd)
#> [1]   4 175
```

You do not need to reshape anything. The scoring functions read the
template’s column names and pull out what they need.

## Scoring delay discounting

[`score_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/score_dd.md)
returns one row per respondent with the discount rate and an attention
flag:

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

- `kval` is the discount rate from Mazur’s hyperbola, mapped from the
  converged indifference delay.
- `ed50` is the effective delay 50%, `1 / k`: the delay at which the
  reward loses half its value. Because it sits on the natural time
  scale, `ed50` is often easier to talk about than `k`.
- `index` records the item at which the adjusting sequence settled.
- `attentionflag` is `"Yes"` when the respondent failed an attention
  check. Those respondents get an `NA` rate, because a failed check
  means the choices cannot be trusted. The third and fourth respondents
  here are flagged for that reason.

### Response times

[`calc_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_dd.md)
returns the full trial-by-trial record, joining each choice to its
response-timing data:

``` r

dd_full <- calc_dd(five.fivetrial_dd)
names(dd_full)
#>  [1] "ResponseId"    "index"         "q"             "firstclick"   
#>  [5] "lastclick"     "pagesubmit"    "totalclicks"   "response"     
#>  [9] "attentionflag" "kval"          "ed50"
head(dd_full, 6)
#> # A tibble: 6 × 11
#>   ResponseId index     q firstclick lastclick pagesubmit totalclicks response
#>        <int> <chr> <dbl>      <dbl>     <dbl>      <dbl>       <dbl> <chr>   
#> 1          1 I16       1       1.76      1.76       3.34           1 ll      
#> 2          1 I24       2       7.73      7.73       8.46           1 ss      
#> 3          1 I20       3       1.56      1.56       3.38           1 ll      
#> 4          1 I22       4       2.33      3.95       4.50           2 ss      
#> 5          1 I21       5       3.16      3.16       3.73           1 ss      
#> 6          2 I16       1       3.78      3.78       4.35           1 ss      
#> # ℹ 3 more variables: attentionflag <chr>, kval <dbl>, ed50 <dbl>
```

Alongside the scored `kval` and `ed50`, you get the trial order (`q`),
the first-click, last-click, and page-submit times, and the click count
for each trial. Response latencies are a useful secondary check on
engagement: implausibly fast submissions often accompany careless
responding.

If you only want the raw choices or the timing on their own,
[`ans_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/ans_dd.md)
and
[`timing_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/timing_dd.md)
return those pieces separately.

## Scoring probability discounting

The probability variant works the same way.
[`score_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/score_pd.md)
reads the probability template and returns a discount rate for odds
against winning:

``` r

data(five.fivetrial_pd)
score_pd(five.fivetrial_pd)
#> # A tibble: 4 × 7
#>   ResponseId index    response  hval attentionflag etheta50  ep50
#>        <int> <chr>    <chr>    <dbl> <chr>            <dbl> <dbl>
#> 1          1 I9       lu        7.44 No              0.134   88.1
#> 2          2 AttendSS lu       99    No              0.0101  99  
#> 3          3 I15      sc        1.60 No              0.624   61.6
#> 4          4 I9       lu        7.44 No              0.134   88.1
```

Here `hval` is the probability-discounting rate, `etheta50` is the
odds-against value at the point of indifference (`1 / h`), and `ep50` is
the corresponding effective probability, `100 / (etheta50 + 1)`,
expressed as a percentage. As with delay discounting,
[`calc_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_pd.md)
adds the trial-level timing, and
[`ans_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/ans_pd.md)
and
[`timing_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/timing_pd.md)
return the components.

## References

- Koffarnus, M. N., & Bickel, W. K. (2014). A 5-trial adjusting delay
  discounting task: Accurate discount rates in less than one minute.
  *Experimental and Clinical Psychopharmacology, 22*(3), 222–228.
  <https://doi.org/10.1037/a0035973>
- Koffarnus, M. N., Kaplan, B. A., & Stein, J. S. (2017). *User guide
  for Qualtrics minute discounting template.*
  <https://doi.org/10.13140/RG.2.2.26495.79527>
- Koffarnus, M. N., Rzeszutek, M. J., & Kaplan, B. A. (2021).
  *Additional discounting rates in less than one minute: Task variants
  for probability and a wider range of delays.*
  <https://doi.org/10.13140/RG.2.2.31281.92000>
