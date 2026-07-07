## Test environments
* locally: Manjaro Linux (kernel 6.6), R 4.6.0
* remotely: Windows Server 2022 x64, R 4.6.0 (win-builder / check_win_release)
* remotely: Ubuntu-latest (GitHub Actions CI)

## R CMD check results

### Version 2.1.0
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

**win-builder (Windows Server 2022, R 4.6.0 release):** Status OK — no errors, warnings, or notes.

**Local (Manjaro Linux, R 4.6.0):** 0 errors, 0 warnings, 1 note:

* `checking compilation flags used`: non-portable flags (`-Werror=format-security`,
  `-march=x86-64`, etc.) originate from the local system compiler configuration,
  not the package. Not present on win-builder or CI.

`cmdstanr` is a suggested package not on CRAN. Its availability via
`Additional_repositories: https://stan-dev.r-universe.dev` was confirmed
by the CRAN incoming check. All cmdstanr-dependent tests are guarded with
`skip_if_not_installed("cmdstanr")` and a check that CmdStan is installed,
so the package checks and tests pass cleanly without CmdStan present.

### Version 2.0.0
0 errors ✔ | 0 warnings ✔ | 2 notes ✖

❯ checking installed package size ... NOTE
    installed size is 261.6Mb
    sub-directories of 1Mb or more:
      libs  260.8Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.


### Version 1.1.0
0 errors ✔ | 0 warnings ✔ | 2 notes ✖

* checking installed package size ... NOTE
    installed size is 10.5Mb
    sub-directories of 1Mb or more:
      libs   9.6Mb

  - Compiled code bmgarch.so is relatively large

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
  
  - GNU make is a build time requirement for Ubuntu
  
## Comments from CRAN admins:
* 06.12.21 and 10.12 Ligges:
  > Pls fix the 404s and 301s from the URL checks.
  - Fixed 404 and 301 errors in README and DESCRIPTION


### Version 1.0.1
* This is a resubmission, initiated by CRAN package check on 2021-06-12

0 Errors or Warnings
2 Notes:

* checking installed package size ... NOTE
    installed size is  9.4Mb
    sub-directories of 1Mb or more:
      libs   8.7Mb

  - Note that compiled code bmgarch.so is 11.5 Mb large

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

  - Note that GNU make is a build time requirement for Ubuntu

## Comments from CRAN admins:
### Version 1.0.0
* 10.09.20 Gregor Seyer: 
> Please replace \dontrun with \donttest
  
  - Done

* 28.08.20 Uwe Ligges: 
> Please single quote software names such as 'rstan' in the Description field.

	- Done


* 26.08.20 The first submission was rejected by Swetlana Herbrandt:

> examples are wrapped in \dontrun{}, hence nothing getst 
> tested. Please unwrap the examples if that is feasible and if they can
> be executed in < 5 sec for each Rd file or create additionally small toy
> examples.

> Alternatively, you can write some tests (e.g. by using testthat). The
> overall execution time should not exceed 10 minutes.

	- Fix: Examples can't be run in less than 5 sec. Instead, we added tests for
	all user-exposed functions for testthat. Overall execution time with 2 cores
	is approx 3.5 mins. Additional examples to illustrate package functions have
	been added as `dontrun{}`

> Please ensure that you do not use more than 2 cores in your examples or
> tests.

	- Done

> If there are references describing the (theoretical background of GARCH
> Models) methods in your package, please add these in the Description
> field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.c

	- Added corresponding information
