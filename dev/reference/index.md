# Package index

## Mixed-effects discounting (TMB)

Fit SLT-beta mixed-effects discounting models to bounded
indifference-point data using TMB, and simulate from the same model.

- [`fit_dd_tmb()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_tmb.md)
  : Fit an indifference-point mixed-effects discounting model via TMB
- [`simulate_dd_ip()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_ip.md)
  : Simulate IP-family mixed-effects discounting data
- [`VarCorr(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/VarCorr.beezdiscounting_tmb.md)
  : Random-effect covariance for a TMB discounting model

## Choice-based discounting (TMB)

Fit the structural or descriptive (Young 2018) SS-vs-LL choice model — a
binomial GLMM — and simulate from either model. The descriptive mode
adds a VarCorr method for the random-effect covariance.

- [`fit_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice.md)
  : Fit a trial-level SS-vs-LL choice model (binomial GLMM) via TMB
- [`simulate_dd_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/simulate_dd_choice.md)
  : Simulate trial-level SS-vs-LL choice data
- [`VarCorr(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/VarCorr.beezdiscounting_choice.md)
  : Random-effect (co)variances for a beezdiscounting_choice fit

## Bayesian discounting (brms)

Bayesian mixed-effects discounting via brms/Stan: the indifference-point
and structural-choice fitters with their default priors.

- [`fit_dd_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_brms.md)
  : Fit a Bayesian Mixed-Effects Discounting Model via brms
- [`fit_dd_choice_brms()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd_choice_brms.md)
  : Fit a Bayesian Structural Choice Discounting Model via brms
- [`default_dd_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_priors.md)
  : Default priors for Bayesian (brms) delay-discounting models
- [`default_dd_choice_priors()`](https://brentkaplan.github.io/beezdiscounting/reference/default_dd_choice_priors.md)
  : Default priors for the Bayesian (brms) choice model

## Estimated marginal means & comparisons

Post-hoc estimated marginal means and pairwise comparisons for fitted
mixed-effects discounting models.

- [`get_dd_param_emms()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_param_emms.md)
  :

  Estimated marginal means of the discount rate `k`

- [`get_dd_comparisons()`](https://brentkaplan.github.io/beezdiscounting/reference/get_dd_comparisons.md)
  :

  Factor-level comparisons of the discount rate `k`

## Model inspection (S3 methods)

S3 methods for extracting information from fitted models.

- [`tidy(`*`<beezdiscounting_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_brms.md)
  : Tidy a beezdiscounting_brms model
- [`tidy(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_choice.md)
  : Tidy a beezdiscounting_choice model into a coefficient tibble
- [`tidy(`*`<beezdiscounting_choice_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_choice_brms.md)
  : Tidy a beezdiscounting_choice_brms model
- [`tidy(`*`<beezdiscounting_comparison>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_comparison.md)
  : Tidy a discounting comparison into a flat contrasts frame
- [`tidy(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/tidy.beezdiscounting_tmb.md)
  : Tidy a beezdiscounting_tmb model into a coefficient tibble
- [`glance(`*`<beezdiscounting_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/glance.beezdiscounting_brms.md)
  : Glance at a beezdiscounting_brms model
- [`glance(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/glance.beezdiscounting_choice.md)
  : Glance at a beezdiscounting_choice model
- [`glance(`*`<beezdiscounting_choice_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/glance.beezdiscounting_choice_brms.md)
  : Glance at a beezdiscounting_choice_brms model
- [`glance(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/glance.beezdiscounting_tmb.md)
  : Glance at a beezdiscounting_tmb model
- [`augment(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/augment.beezdiscounting_choice.md)
  : Augment a beezdiscounting_choice model
- [`augment(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/augment.beezdiscounting_tmb.md)
  : Augment a beezdiscounting_tmb model
- [`predict(`*`<beezdiscounting_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/predict.beezdiscounting_brms.md)
  : Predict from a beezdiscounting_brms model
- [`predict(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/predict.beezdiscounting_choice.md)
  : Predict from a structural choice discounting model
- [`predict(`*`<beezdiscounting_choice_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/predict.beezdiscounting_choice_brms.md)
  : Predict P(LL) from a beezdiscounting_choice_brms model
- [`predict(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/predict.beezdiscounting_tmb.md)
  : Predict from a TMB mixed-effects discounting model
- [`confint(`*`<beezdiscounting_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/confint.beezdiscounting_brms.md)
  : Credible intervals for a beezdiscounting_brms model
- [`confint(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/confint.beezdiscounting_choice.md)
  : Confidence intervals for a structural choice discounting model
- [`confint(`*`<beezdiscounting_choice_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/confint.beezdiscounting_choice_brms.md)
  : Credible intervals for a beezdiscounting_choice_brms model
- [`confint(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/confint.beezdiscounting_tmb.md)
  : Confidence intervals for a TMB discounting model
- [`summary(`*`<beezdiscounting_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/summary.beezdiscounting_brms.md)
  : Summarize a beezdiscounting_brms model
- [`summary(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/summary.beezdiscounting_choice.md)
  : Summarize a structural choice discounting fit
- [`summary(`*`<beezdiscounting_choice_brms>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/summary.beezdiscounting_choice_brms.md)
  : Summarize a beezdiscounting_choice_brms model
- [`summary(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/summary.beezdiscounting_tmb.md)
  : Summarize a TMB mixed-effects discounting fit
- [`coef(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/coef.beezdiscounting_choice.md)
  : Extract coefficients from a structural choice model
- [`coef(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/coef.beezdiscounting_tmb.md)
  : Extract coefficients from a TMB discounting model
- [`fixef(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/fixef.beezdiscounting_choice.md)
  : Extract fixed effects from a structural choice model
- [`fixef(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/fixef.beezdiscounting_tmb.md)
  : Extract fixed effects from a TMB discounting model
- [`ranef(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/ranef.beezdiscounting_choice.md)
  : Extract subject-level random effects from a choice model
- [`ranef(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/ranef.beezdiscounting_tmb.md)
  : Extract subject-level random effects from a TMB discounting model
- [`fitted(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/fitted.beezdiscounting_choice.md)
  : Fitted values for a beezdiscounting_choice fit
- [`fitted(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/fitted.beezdiscounting_tmb.md)
  : Fitted values for a beezdiscounting_tmb fit
- [`residuals(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/residuals.beezdiscounting_choice.md)
  : Residuals for a beezdiscounting_choice fit
- [`residuals(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/residuals.beezdiscounting_tmb.md)
  : Residuals for a beezdiscounting_tmb fit
- [`print(`*`<beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/print.beezdiscounting_choice.md)
  : Print a structural choice discounting fit
- [`print(`*`<beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/print.beezdiscounting_tmb.md)
  : Print a TMB mixed-effects discounting fit
- [`print(`*`<summary.beezdiscounting_choice>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/print.summary.beezdiscounting_choice.md)
  : Print a structural choice discounting model summary
- [`print(`*`<summary.beezdiscounting_tmb>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/print.summary.beezdiscounting_tmb.md)
  : Print a TMB discounting model summary

## Delay discounting (NLS / scoring)

Fit hyperbolic/exponential discount functions via nonlinear regression
and score delay-discounting tasks.

- [`fit_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/fit_dd.md)
  : Fit Delay-Discounting Model
- [`results_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/results_dd.md)
  : Extract Results from Delay-Discounting Model
- [`plot_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/plot_dd.md)
  : Plot Delay-Discounting Model
- [`calc_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_dd.md)
  : Calculate scores, answers, and timing for 5.5 trial delay
  discounting from Qualtrics template
- [`score_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/score_dd.md)
  : Score 5.5 trial delay discounting from Qualtrics template
- [`ans_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/ans_dd.md)
  : Converts answers from 5.5 trial delay discounting from Qualtrics
  template
- [`timing_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/timing_dd.md)
  : Extract timing metrics from 5.5 trial delay discounting from
  Qualtrics template

## Probability discounting

Scoring functions for probability-discounting tasks.

- [`calc_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_pd.md)
  : Calculate scores, answers, and timing for 5.5 trial probability
  discounting from Qualtrics template
- [`score_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/score_pd.md)
  : Score 5.5 trial probability discounting from Qualtrics template
- [`ans_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/ans_pd.md)
  : Converts answers from 5.5 trial probability discounting from
  Qualtrics template
- [`timing_pd()`](https://brentkaplan.github.io/beezdiscounting/reference/timing_pd.md)
  : Extract timing metrics from 5.5 trial probability discounting from
  Qualtrics template

## 27-item Monetary Choice Questionnaire (MCQ)

Scoring, lookup, simulation, and reshaping for the 27-item MCQ.

- [`score_mcq27()`](https://brentkaplan.github.io/beezdiscounting/reference/score_mcq27.md)
  : Score 27-item MCQ
- [`summarize_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/summarize_mcq.md)
  : Provide a summary of the results from the MCQ output table.
- [`prop_ss()`](https://brentkaplan.github.io/beezdiscounting/reference/prop_ss.md)
  : Calculate proportion of SIR/SS responses at each k value
- [`get_lookup_table()`](https://brentkaplan.github.io/beezdiscounting/reference/get_lookup_table.md)
  : Get internal lookup table for the 27-item MCQ
- [`mcq27_to_choice()`](https://brentkaplan.github.io/beezdiscounting/reference/mcq27_to_choice.md)
  : Convert 27-item MCQ responses to a trial-level choice frame
- [`generate_data_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/generate_data_mcq.md)
  : Generate fake MCQ data
- [`wide_to_long_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/wide_to_long_mcq.md)
  : Reshape MCQ data wide to long
- [`long_to_wide_mcq()`](https://brentkaplan.github.io/beezdiscounting/reference/long_to_wide_mcq.md)
  : Reshape MCQ data long to wide
- [`wide_to_long_mcq_excel()`](https://brentkaplan.github.io/beezdiscounting/reference/wide_to_long_mcq_excel.md)
  : Reshape MCQ data from wide (as used in the 21- and 27-Item Monetary
  Choice Questionnaire Automated Scorer) to long
- [`long_to_wide_mcq_excel()`](https://brentkaplan.github.io/beezdiscounting/reference/long_to_wide_mcq_excel.md)
  : Reshape MCQ data from long to wide (as used in the 21- and 27-Item
  Monetary Choice Questionnaire Automated Scorer)

## Utilities & metrics

Area under the curve, fit metrics, and systematicity checks.

- [`calc_aucs()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_aucs.md)
  : Calculate Area-Under-the-Curve (AUC) Metrics for Delay Discounting
  Data
- [`calc_r2()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_r2.md)
  : Calculate R-Squared for a Model
- [`calc_conf_int()`](https://brentkaplan.github.io/beezdiscounting/reference/calc_conf_int.md)
  : Calculate Confidence Intervals for a Parameter
- [`check_unsystematic()`](https://brentkaplan.github.io/beezdiscounting/reference/check_unsystematic.md)
  : Check for Unsystematic Data Violations

## Visualization

Plotting methods for scored task outputs.

- [`plot(`*`<prop_ss_output>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/plot.prop_ss_output.md)
  : Plot Proportion of SIR/SS Choices by k Value
- [`plot(`*`<score_mcq27_output>`*`)`](https://brentkaplan.github.io/beezdiscounting/reference/plot.score_mcq27_output.md)
  : Plot MCQ-27 Scores
- [`plot_dd()`](https://brentkaplan.github.io/beezdiscounting/reference/plot_dd.md)
  : Plot Delay-Discounting Model

## Datasets

Built-in datasets for examples and testing.

- [`dd_ip`](https://brentkaplan.github.io/beezdiscounting/reference/dd_ip.md)
  : Delay Discounting Data
- [`mcq27`](https://brentkaplan.github.io/beezdiscounting/reference/mcq27.md)
  : Example 27-item MCQ data
- [`five.fivetrial_dd`](https://brentkaplan.github.io/beezdiscounting/reference/five.fivetrial_dd.md)
  : Example Qualtrics output from the 5.5 trial delay discounting
  template.
- [`five.fivetrial_pd`](https://brentkaplan.github.io/beezdiscounting/reference/five.fivetrial_pd.md)
  : Example Qualtrics output from the 5.5 trial probability discounting
  template.

## Reexports

Generics re-exported from other packages.

- [`reexports`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  [`augment`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  [`tidy`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  [`glance`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  [`fixef`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  [`ranef`](https://brentkaplan.github.io/beezdiscounting/reference/reexports.md)
  : Objects exported from other packages
- [`` `%>%` ``](https://brentkaplan.github.io/beezdiscounting/reference/pipe.md)
  : Pipe operator
