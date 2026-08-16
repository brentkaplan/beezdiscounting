## Submission

This is a feature release (0.4.0) of an existing CRAN package (0.3.2, January
2025). It adds:

* a mixed-effects (hierarchical) indifference-point discounting model fit via
  Template Model Builder ('TMB') under a scale-location-truncated beta
  likelihood that admits indifference points at exactly 0 and 1 (`fit_dd_tmb()`),
  with the accompanying broom / 'emmeans'-style S3 surface, group comparisons,
  and a simulator; the two-parameter Green-Myerson and Rachlin hyperboloids
  and a subject-random curvature exponent;
* trial-level choice models (`fit_dd_choice()`; structural and Young-2018
  descriptive modes);
* an optional Bayesian tier via 'brms' (`fit_dd_brms()`,
  `fit_dd_choice_brms()`); 'brms', 'posterior' and 'loo' are Suggests only and
  every example/test/vignette that needs them is conditional on their presence;
* Monte Carlo power analysis (`power_discounting()`, `find_n_discounting()`);
* scoring of the 21-item Monetary Choice Questionnaire and the 30-item
  Probability Discounting Questionnaire, plus bundled example data;
* ten vignettes (all new since 0.3.2, which shipped none).

Compiled code: the 'TMB' C++ templates in `src/` changed since 0.3.2 (new
model families and random-effect blocks).

## Test environments

* local macOS (R 4.5.2)
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest,
  windows-latest (R-CMD-check workflow)
* win-builder (R-release and R-devel) and R-hub v2 (linux, macos-arm64,
  windows, clang-asan, valgrind -- the package ships TMB C++ templates)

## R CMD check results

Locally (macOS/arm64, R 4.5.x): 0 errors | 1 warning | 1 note -- the single
local warning is environmental (R's own header under a very recent Apple
'clang'), explained below.

win-builder R-release (R 4.6.1 ucrt, Windows Server 2022, 2026-08-16):
0 errors | 0 warnings | 1 note (the pkgdown-site URL note below).
win-builder R-devel (R Under development 2026-08-15 r90413 ucrt, 2026-08-16):
0 errors | 0 warnings | 1 note (the pkgdown-site URL note below).
R-hub v2 (2026-08-16; R-devel on linux/macos-arm64/windows, plus clang-asan
and valgrind): all five legs OK -- 0 errors, no sanitizer findings; valgrind
`ERROR SUMMARY: 0 errors`, `definitely lost: 0 bytes`.

* **NOTE — possibly invalid URL** `https://brentkaplan.github.io/beezdiscounting/`:
  this is the package's 'pkgdown' documentation site, which is (re)deployed by the
  release that accompanies this submission; the URL resolves once the site is
  published.

* **Install size (~13 Mb):** the package links to 'TMB' / 'RcppEigen' and ships a
  compiled shared library; the `libs/` directory accounts for the size, as is
  typical for a compiled model package.

* **Local-only install WARNING** `-Wfixed-enum-extension`: this originates in R's
  own header (`R_ext/Boolean.h`) under a very recent Apple 'clang', not in package
  code, and does not appear on CRAN's build machines.

* The DOIs in the README that automated checkers report as "403 Forbidden"
  (e.g. ResearchGate-hosted DOIs) are valid registered DOIs that block automated
  HEAD requests; they resolve in a browser.

## Method references

References for the methods are given in the Description field, including the new
mixed-effects scale-location-truncated beta discounting model
(Kim, Kaplan, Koffarnus, and Franck, 2025; <arXiv:2509.13167>).

## revdepcheck results

There are currently no downstream dependencies for this package.
