## Submission

This is a new submission.

## Test environments

* Local: Ubuntu 24.04 R 4.3.3 + Fedora 44 R 4.6.2, with the suggested packages installed
  (`R CMD check --as-cran`).
* GitHub Actions: ubuntu-latest, R release.
* winbuilder: all 3 options.
* The C++ code was also run through the whole test suite built with
  AddressSanitizer and UndefinedBehaviorSanitizer (gcc 13), without any error.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: Emik, Terrill. These are the
  names of the authors of the method implemented by the package (Emik and
  Terrill 1949, cited with its DOI), and are spelled correctly.
