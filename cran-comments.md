## Submission

This is a feature release (0.4.0) of an existing CRAN package. It adds a
mixed-effects (hierarchical) indifference-point discounting model fit via
Template Model Builder ('TMB'), under a scale-location-truncated beta likelihood
that admits indifference points at exactly 0 and 1, with the accompanying broom /
'emmeans'-style S3 surface, a simulator, and a vignette.

## Test environments

* local macOS (R 4.5.2)
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest,
  windows-latest (R-CMD-check workflow)

## R CMD check results

0 errors | 0 warnings | 1 note (locally an additional environmental install
warning, explained below).

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
