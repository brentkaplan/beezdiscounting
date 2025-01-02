# beezdiscounting 0.3.2

### New Features

- **`fit_dd()`**:
  - Introduced a new function to fit delay-discounting models using specified equations (`"mazur"`/`"hyperbolic"` or `"exponential"`) and methods (`"pooled"`, `"mean"`, or `"two stage"`).
  - Supports flexible data handling for aggregated and participant-specific modeling.
  - Returns an object of class `"fit_dd"` containing the fitted models, input data, and method details.

- **`plot_dd()`**:
  - Added a function to visualize fitted delay-discounting models.
  - Automatically adapts to different fitting methods, including aggregated and individual models.
  - Provides customizable axis labels, title, and optional log-transformed x-axis for improved visualization of delay scales.

- **`results_dd()`**:
  - New utility to extract model parameter estimates, confidence intervals, and fit statistics from a `"fit_dd"` object.
  - Supports both aggregated and participant-specific models.
  - Outputs a tidy tibble with columns for terms, estimates, standard errors,
    t-statistics, p-values, R2, three different AUC metrics, and confidence bounds.

- **`check_unsystematic()`**:
  - New utility function to check delay-discounting datasets for unsystematic
    data patterns according to Johnson & Bickel's (2008) two criteria.

- **`calc_aucs()`**:
  - New utility function to calculate three different area under the curve
    (AUC) metrics for delay-discounting data according to Borges et al. (2016).

### Improvements

- Confidence intervals are now computed using the `calc_conf_int()` function, ensuring accurate estimation based on model degrees of freedom.
- R2 values are calculated consistently using the `calc_r2()` function, providing reliable fit metrics for all models.

### Enhancements

- The package now supports robust delay-discounting workflows, from unsystematic
  identification (`check_unsystematic`), model fitting (`fit_dd`), to visualization (`plot_dd`), to result extraction
  (`results_dd`).
- Improved compatibility with delay-discounting datasets that require participant-level or aggregated modeling approaches.


# beezdiscounting 0.3.1

## Minor fix

* Correctly names output columns from `calc_pd()` and `score_pd()`. `ep50` changed to `etheta50` and corrected calculation of `ep50`.

# beezdiscounting 0.3.0

## New features

* Add functions for scoring 5.5 trial probability discounting task (from the Qualtrics template) including: `calc_pd()`
(and `score_pd()`, `timing_pd()`, and `ans_pd`).

## Minor fix

* Subsetting issue is fixed in `score_dd()` that would unintentionally drop all rows if both conditions were `FALSE`.

## Other changes

* Rename example data from `five.fivetrial` to `five.fivetrial_dd` for delay discounting.

* Add example data `five.fivetrial_pd` for probability discounting.

# beezdiscounting 0.2.0

## New features

* `score_mcq27()` properly supports arguments: `impute_method`, `random`, `return_data`, and `verbose`.
See documentation and the `README` for explanations.

* `generate_data_mcq()` can generate fake MCQ data, including `seed` and `prop_na` arguments for
reproducibility and specifying proportion of `NA`s.

* `long_to_wide*` and `wide_to_long*` are helper functions to reshape data from/to different formats.

## Minor fix

* When no imputation is specified and `NA`s exist in the data, `score_mcq27()` returns `NA`s for the scoring
instead of 1.

# beezdiscounting 0.1.0

* Initial release with basic scoring of 27-item Monetary Choice Questionnaire and 5.5 trial delay discounting task from the Qualtrics template.

* Added a `NEWS.md` file to track changes to the package.
