## R CMD check results

- running `devtools::check(remote = TRUE, manual = TRUE)` gives no errors, warnings, and 1 note:
```
❯ checking CRAN incoming feasibility ... [4s/58s] NOTE
  Maintainer: ‘Tomasz Woźniak <wozniak.tom@pm.me>’
  
  New submission
```
We are happy about it.

- the package passes all the checks on GH action from the standard setup `R-CMD-check.yaml`

- the package passes all procedures following `usethis::use_release_issue()`

- comments from CRAN
```
Thanks,

Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'C++'
Please note that package names are case sensitive.
For more details:
<https://contributor.r-project.org/cran-cookbook/description_issues.html#formatting-software-names>

Please omit the redundant "for R Packages" at the end/start of your
title and description.

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")
For more details:
<https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

Please fix and resubmit.

Best,
Benjamin Altmann
```
We implemented all changes except for *If there are references...*. In an email exchange, we agreed this is OK. Thanks!