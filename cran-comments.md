## Test environments
* local ubuntu 18.10 install, R 4.0.2
* ubuntu 16.04.6 LTS (on travis-ci), R Under development (unstable) (2020-08-23 r79071)
* MacOS-latest and Windows-latest (on github acttions-ci)
* win_devel

  
## R CMD check restults

0 Errors or Warnnings
2 Notes:

* checking installed package size ... NOTE
    installed size is  8.8Mb
    sub-directories of 1Mb or more:
      libs   8.1Mb

  Compiled code bmgarch.so is 8.1 Mb large

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

  GNU make is a build time requirement for Ubuntu

* This is a resubmission of the first version of the package to CRAN

## Comments from CRAN admins:
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
