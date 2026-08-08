# Internal MCQ version registry

Internal MCQ version registry

## Usage

``` r
.mcq_registry(items)
```

## Arguments

- items:

  Number of MCQ items (27 or 21)

## Value

List with the version's design table (in ladder order), items,
per-magnitude item count, overall ladder edge k, and k-rank labels.

## Details

Edge conventions mirror the Kaplan et al. (2014) Excel scorers: overall
ladders append 0.25 (27-item) / 0.1333 (21-item) past the steepest item;
magnitude ladders repeat their last kindiff.
