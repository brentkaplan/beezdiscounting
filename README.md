# Behavioral Economic (be) Easy (ez) Discounting
An R package containing commonly used functions for analyzing behavioral economic discounting data.

The package supports scoring of the 27-Item Monetary Choice
Questionnaire (see Kaplan et al., 2016; <doi:10.1007/s40614-016-0070-9>) and scoring of the 
minute discounting task (see Koffarnus & Bickel, 2014; <doi:10.1037/a0035973>) using the 
Qualtrics 5-trial discounting template (see the Qualtrics Minute Discounting User Guide; 
<doi:10.13140/RG.2.2.26495.79527>), which is also available as a .qsf file in this package.

## Note About Use
Currently in development.

## Installing beezdemand

### GitHub Release

To install a stable release directly from
[GitHub](https://github.com/brentkaplan/beezdiscounting), first install and
load the `devtools` package. Then, use `install_github` to install the
package and associated vignette. You *don’t* need to download anything
directly from [GitHub](https://github.com/brentkaplan/beezdiscounting), as
you should use the following instructions:

``` r
install.packages("devtools")

devtools::install_github("brentkaplan/beezdiscounting", build_vignettes = TRUE)

library(beezdiscounting)
```

### Questions, Suggestions, and Contributions
---------------------------------------------

Have a question? Have a suggestion for a feature? Would you like to contribute? Email me at <bkaplan.ku@gmail.com>.

### License
-----------

GPL (>= 2)
