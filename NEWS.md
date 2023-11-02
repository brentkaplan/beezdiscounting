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
