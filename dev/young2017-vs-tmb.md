# Young (2017) vs. beezdiscounting's TMB indifference-point tier

A point-by-point comparison of the multilevel approach Young recommends for
indifference-point discounting data against what `fit_dd_tmb()` actually does.

- **Paper:** Young, M. E. (2017). Discounting: A practical guide to multilevel
  analysis of indifference data. *Journal of the Experimental Analysis of
  Behavior, 108*(1), 97–112. doi:[10.1002/jeab.265](https://doi.org/10.1002/jeab.265).
  Companion (choice data, *not* this comparison): Young (2018),
  doi:[10.1002/jeab.316](https://doi.org/10.1002/jeab.316).
- **Our side:** `R/dd-tmb.R` (wrapper) and `src/MixedDiscounting.h` (the C++ TMB
  template — the authoritative model form). Line anchors below are as of this
  writing.
- **One-line takeaway:** `fit_dd_tmb()` *is* Young's one-stage multilevel
  recommendation, implemented on TMB (a Laplace approximation with exact automatic
  differentiation), and extended with a bounded-support error family (SLT-beta)
  that addresses the heteroscedasticity and boundary problems Young flags but
  leaves unsolved.

## Comparison table

| # | Axis | Young (2017) | `fit_dd_tmb()` | Verdict |
|---|------|--------------|----------------|---------|
| 1 | **Discount function(s)** | Hyperbolic/Mazur (Eq. 1), and the two-parameter hyperboloids Green-Myerson (Eq. 3) and Rachlin (Eq. 4). Frames the methods as applying to "any discounting function" (p. 97) but demonstrates these three; exponential is not among the fitted functions. | `mazur`, `green-myerson`, `rachlin` — **plus `exponential`** (`eqn_type` 0/1/2/3, `R/dd-tmb.R:161`; math `MixedDiscounting.h:11–14, 171–189`). | **Match + extend** (we add exponential) |
| 2 | **Parameterization & link** | Estimates **`logk`, not `k`** ("logk value will be estimated rather than `k`", p. 98–99; his Eq. 2 = `Amount / (1 + exp(logk)·Delay)`). Response is the indifference **amount** (value scale). | Linear predictor `log k = X·beta_k`; `k = exp(·)` (`MixedDiscounting.h:18–21, 144`). Response is the **proportion** `y ∈ [0,1]` (amount coerced via `ll`/`response_scale`). Identity link from the discount function to the mean. | **Match** (same log-rate idea; we normalize to a proportion) |
| 3 | **Random-effects structure** | Random intercept on `logk` across subjects (1-param). For the hyperboloids, **correlated** random effects on `(logk, s)`; reports the RE correlation (Green-Myerson r = −.87, Rachlin r = −.39, p. 108) and flags GM over-parameterization. | `k ~ 1`: random intercept on `log k` (`MixedDiscounting.h:109–116`). `k + phi ~ 1` / `k + s ~ 1`: a 2-D random intercept via a 2×2 Cholesky (correlated under `pdSymm`, independent under `pdDiag`; `rho = tanh(cor_re)`) (`MixedDiscounting.h:117–134`; `R/dd-tmb.R:1144–1147, 1265–1287`). | **Match** for `(log k, s)`; we **add** a random-precision option `(log k, phi)` |
| 4 | **Estimation / inference** | One-stage MLE via **R `nlme`** (Lindstrom–Bates alternating algorithm; first-order conditional linearization of the nonlinear model), seeded by `nlsList` individual fits (p. 100). Wald SEs from `nlme`. | One-stage MLE via **TMB**: `MakeADFun(..., random=)` integrates REs by **Laplace approximation** with exact autodiff; `nlminb` (default) / `optim(L-BFGS-B)`; **3-set multi-start** kept on a sanity guard (`R/dd-tmb.R:684–816`); SEs via `sdreport` with a `pdHess` gate (`R/dd-tmb.R:867–895`). | **Same goal, different engine** (exact-AD Laplace vs. nlme's first-order linearization) |
| 5 | **Error / likelihood** | **Normal/Gaussian** residuals on the indifference value. Explicitly names "normality of the residuals and homogeneity of variance" as the live assumptions, notes observed **heteroscedasticity**, and says addressing it "can be complex" (p. 110–111). | **SLT-beta by default** (bounded-support (0,1) beta; mean = discount function, precision `phi`; variance ≈ `mu(1−mu)/(1+phi)` for the near-identity SLT transform, hence mean-dependent → heteroscedasticity handled automatically; `MixedDiscounting.h:196–228`). `gaussian` offered as the Young-equivalent option (`:230–237`). | **Primary divergence / our value-add** |
| 6 | **Points at exactly 0 / 1** | Not addressed. Normal errors on the value scale admit impossible predictions and treat boundary piles as ordinary points. | SLT-beta admits exact 0 and 1 via a **scale-location-truncation** transform (`y/s_slt + l`) and a truncation normalizer `Z` (`MixedDiscounting.h:27–34, 219–227`). | **We solve a problem Young doesn't raise** (Kim et al.) |
| 7 | **Covariates / group effects** | Group/condition differences as **fixed effects on `logk`** inside the one model. Strongly recommends **effect coding** and **centering** with interactions (p. 110–111); treat graded vars (day/block/session) as continuous single-df. | Between-subject `factors` (+ pairwise `factor_interaction`) and `continuous_covariates` enter `X` on `log k`; group contrasts via `emmeans` (`get_dd_param_emms()` / `get_dd_comparisons()`). Design uses `model.matrix` with the session's contrasts (treatment by default), no auto-centering (`R/dd-tmb.R:106–140`). | **Match in structure; we lag on coding guidance** (see follow-ups) |
| 8 | **Model selection & diagnostics** | One fit-statistic over all data; **AIC/BIC** preferred (p. 102), **AICc** in simulations. "Always examine a plot of the model superimposed on the data"; residuals-vs-fitted for heteroscedasticity. Convergence problems discussed at length. | `glance()` → logLik/**AIC/BIC** (marginal; counts fixed params only, REs integrated out, `R/dd-tmb.R:1444–1448`). `equation`/`family` are user-chosen, not auto-selected. `plot()` / `plot_qq()` diagnostics; multi-start + `pdHess` gate + blow-up/phi guards for the convergence pitfalls Young describes. | **Match** (we expose AIC/BIC, not AICc) |
| 9 | **One-stage vs. two-stage** | The paper's **central thesis**: reject fit-per-subject-then-average (ignores differential precision, mishandles missing data, invites exclusion-criterion researcher d.f., 2-param per-subject fits are unstable). One-stage MLM pools strength; **shrinkage** toward the group; "operates like a Bayesian analysis" with the group as a prior (p. 97–102). | `fit_dd_tmb()` is exactly that one-stage hierarchical fit — group (fixed) and subject (random) estimated jointly from one marginal likelihood; subject `k` are shrunken BLUP-style reconstructions (`R/dd-tmb.R:1027–1075`). | **Full alignment** (this is the shared philosophy) |
| 10 | **Practical recommendations** | Estimate `logk`; seed from a simpler fit; plot model over data; effect-code + center; prefer **Rachlin over Green-Myerson** (GM over-parameterized, p. 108–109); expect implausible per-subject `s`>1 that MLM tames; multilevel as default. | We do: `logk`, data-driven multi-start, `plot()`, multilevel-by-default, bounded `s` (optimizer bounds on population `s`; smooth clamp on subject-level `s`). We don't (yet): effect coding/centering, an in-docs Rachlin-over-GM steer, AICc. | **Mostly followed; a few doc gaps** |

## Where we align with Young

The backbone is identical, and deliberately so. `fit_dd_tmb()` is a **one-stage
nonlinear multilevel model** that estimates the group and the subjects together
from a single likelihood, precisely the approach Young argues for throughout the
paper over the two-stage fit-then-aggregate workflow. Both estimate the discount
rate on the **log scale** (`log k`), for the same reasons Young gives (convergence
and the skew of `k`; his Eq. 2 is our linear predictor). Both place a **random
intercept on `log k`** by default, and both extend to **correlated random effects
on the two hyperboloid parameters** `(log k, s)`. Young even reports the
`log k`–`s` random-effect correlation that our `k + s ~ 1` (`pdSymm`) mode
estimates as a free parameter. The shrinkage Young praises (extreme or sparse
subjects pulled toward the group, with no consistency-criterion culling) is a
direct consequence of our random effects. Young's practical-fitting advice also
maps over cleanly: estimate `log k`, seed the optimizer sensibly, plot the model
over the data, and expect two-parameter convergence trouble that the multilevel
structure helps stabilize.

## Where we go beyond Young

**The error distribution is where the two approaches most substantively diverge.**
Young fits **normal** residuals on the indifference value. He is explicit that the
two key assumptions are normality and homogeneity of variance, observes that
discounting residuals are heteroscedastic, and concedes that the fix "can be
complex" (p. 110–111) without pursuing it. Our **SLT-beta** family provides that
treatment. The beta distribution has bounded support and a mean-dependent variance
(approximately `mu(1−mu)/(1+phi)`, since the SLT scale-location-truncation shifts it
slightly), so the non-constant variance Young flags is modeled rather than assumed
away. Because `mu` is the discount function, the residual spread changes with delay
on its own (largest near `mu = 0.5`, shrinking toward the bounds), with no separate
variance function to specify; the dispersion `phi` stays constant across delays
(per-subject under `k + phi ~ 1`). This is the gap Kim, Koffarnus & Franck (2024,
*Thinking Inside the Bounds*, arXiv:2404.18000) set out to close. Their paper names
the normal, constant-variance model as the baseline and beta regression as the
remedy for delay-dependent variance and out-of-bounds predictions.

**Boundary observations at exactly 0 and 1.** Young's normal model treats the
observations at the bounds as ordinary points (and can predict impossible values).
Plain beta regression is undefined at exactly 0 or 1, where the log-likelihood
diverges. The SLT-beta **scale-location-truncation** trick (Kim, Kaplan, Koffarnus
& Franck, 2025, arXiv:2509.13167) admits exact 0 and 1 without discarding or nudging
them, a problem Young never raises but that real titration data force.

**A random-precision option.** `k + phi ~ 1` lets each subject carry their own
SLT-beta precision (a per-subject noise level). Young has no analog, given that a
normal model has a single residual variance and no precision parameter to
randomize.

**A different estimation engine.** `nlme` marginalizes the random effects by the
Lindstrom–Bates first-order linearization; TMB uses a true **Laplace approximation
with exact automatic-differentiation gradients and Hessians**. Both target the same
one-stage MLE; TMB's Laplace approximation is a higher-order method than nlme's
first-order linearization, computed with exact-AD gradients, and is run here with
multi-start and explicit degeneracy guards (the phi floor, the log-k blow-up bound,
and a smooth clamp on subject-level `s`).

**A Bayesian tier.** Young notes that his MLE multilevel fit "operates like a
Bayesian analysis." We also provide an actual Bayesian tier (brms/Stan) for users
who want a fully Bayesian treatment.

## Where we lag Young (actionable)

1. **Effect coding and centering.** Young strongly recommends effect coding and
   centering covariates, especially with interactions, to cut nonessential
   collinearity and make coefficients interpretable (p. 110–111). Our design calls
   `model.matrix()` with the session's contrasts (treatment by default) and does not
   center. In practice the `emmeans` layer makes the *marginal* means and contrasts
   coding-invariant, so the conclusions are safe; the raw `beta_k` coefficients and
   any interaction terms, however, inherit the reference-level coding. *Follow-up:* a
   docs note (emmeans contrasts are coding-invariant, with centering guidance), or
   optional support for `contr.sum`.
2. **Rachlin-over-Green-Myerson steer.** Young concludes that researchers using a
   hyperboloid "are advised to use multilevel modeling of the Rachlin version"
   because Green-Myerson is over-parameterized (RE `log k`–`s` r = −.87 vs. −.39,
   and a better AIC, p. 108–109). We expose both equally with no guidance.
   *Follow-up:* a sentence in `?fit_dd_tmb` and the TMB vignette. (Our `pdSymm`
   two-RE mode surfaces that same `log k`–`s` correlation, so a user can inspect
   the Green-Myerson problem for themselves.)
3. **AICc.** Young uses small-sample-corrected AIC; we report AIC/BIC only. Minor;
   add AICc to `glance()` if desired.

## Doc/code follow-ups (from this memo)

Applied in this change:

- Added `@references Young (2017)` (plus the Kim et al. SLT-beta sources) to
  `fit_dd_tmb()`'s roxygen.
- Added a "Relationship to Young (2017)" section to `vignette("tmb-mixed-effects")`
  (same one-stage model, bounded-error likelihood beyond it).
- Added a factor-coding note to the TMB vignette's practical notes (emmeans
  contrasts are coding-invariant; effect coding and centering aid interpretation).

Still open:

- A Rachlin-over-Green-Myerson steer in `?fit_dd_tmb` (Young advises the Rachlin
  hyperboloid; Green-Myerson is over-parameterized).
- AICc in `glance()` (Young uses the small-sample correction; we report AIC/BIC).
- The TMB-vignette intro cites Young (2017) for "this style of multilevel analysis."
  That wording is accurate and can stay; the model genuinely implements Young's
  one-stage NLME with a different error family, so the citation is, if anything,
  conservative.

## References

- Young, M. E. (2017). Discounting: A practical guide to multilevel analysis of
  indifference data. *J Exp Anal Behav, 108*(1), 97–112.
  doi:[10.1002/jeab.265](https://doi.org/10.1002/jeab.265).
- Kim, M., Koffarnus, M. N., & Franck, C. T. (2024). *Thinking Inside the Bounds:
  Improved Error Distributions for Indifference Point Data Analysis and Simulation
  via Beta Regression Using Common Discounting Functions.* arXiv:2404.18000.
- Kim, M., Kaplan, B. A., Koffarnus, M. N., & Franck, C. T. (2025).
  *Scale-Location-Truncated Beta Regression: Expanding Beta Regression to
  Accommodate 0 and 1.* arXiv:2509.13167.

*Paper metadata via PubMed (PMID 28699271); full text supplied by the author.*
