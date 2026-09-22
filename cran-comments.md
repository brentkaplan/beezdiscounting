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
* ten vignettes (all new since 0.3.2, which had none).

Compiled code: the 'TMB' C++ templates in `src/` changed since 0.3.2 (new
model families and random-effect blocks).

## Test environments

* local macOS 26 / arm64 (R 4.5.2), `R CMD check --as-cran --run-donttest`
* win-builder R-release (R 4.6.1 ucrt, Windows Server 2022, 2026-09-22)
* win-builder R-devel (R Under development 2026-09-21 r90579 ucrt, 2026-09-22)
* R-hub v2 (2026-09-22): linux (R-devel), windows (R-devel), clang-asan,
  valgrind -- the package contains TMB C++ templates
* GitHub Actions (2026-09-22): ubuntu-latest (R-release, R-devel, oldrel-1),
  windows-latest, macOS 15 and macOS 26 (R-release)
* GitHub Actions, CRAN-like configuration (2026-09-22; `--as-cran`, `NOT_CRAN`
  unset): ubuntu-latest, macOS 15, macOS 26, windows-latest (R-release)

## R CMD check results

win-builder R-release: 0 errors | 0 warnings | 1 note.
win-builder R-devel: 0 errors | 0 warnings | 1 note.
Both notes are the incoming-feasibility note giving the maintainer name and
"GNU make is a SystemRequirements" (see below).

R-hub v2 (2026-09-22): linux, windows, clang-asan and valgrind all OK -- 0
errors, no sanitizer findings; valgrind `ERROR SUMMARY: 0 errors`,
`definitely lost: 0 bytes`. The macos-arm64 leg could not be run: its package
repository currently serves no macOS arm64 binaries for R-devel, so dependency
installation fails before the check starts. macOS is covered by the four macOS
checks on GitHub Actions above (R-release on macOS 15 and macOS 26, one of them
under the CRAN-like configuration).

GitHub Actions: all legs OK in both workflows.

Locally (macOS/arm64, R 4.5.2): 0 errors | 1 warning | 1 note. The warning is
environmental (R's own header under a very recent Apple 'clang') and the note
is one `\donttest` example just over five seconds; both are explained below.

* **NOTE, GNU make is a SystemRequirements:** 'TMB' packages require GNU make
  to build the C++ templates; the dependency is declared in `SystemRequirements`.

* **Install size (~13 Mb)** where it is reported: the package links to 'TMB' /
  'RcppEigen' and installs a compiled shared library; the `libs/` directory
  accounts for the size, as is typical for a compiled model package.

* **Local-only install WARNING** `-Wfixed-enum-extension`: this originates in R's
  own header (`R_ext/Boolean.h`) under a very recent Apple 'clang', not in package
  code, and does not appear on win-builder, R-hub or CRAN's build machines.

* **Local-only example-timing NOTE:** `find_n_discounting` runs 5.3 s on the
  submitter's machine. The example is inside `\donttest{}` and was measured with
  `--run-donttest`; win-builder reports the examples at 14 s in total.

* The DOIs in the README that automated checkers report as "403 Forbidden"
  (e.g. ResearchGate-hosted DOIs) are valid registered DOIs that block automated
  HEAD requests; they resolve in a browser.

## Method references

References for the methods are given in the Description field, including the new
mixed-effects scale-location-truncated beta discounting model
(Kim, Kaplan, Koffarnus, and Franck, 2025; <arXiv:2509.13167>).

## revdepcheck results

There are currently no downstream dependencies for this package.
