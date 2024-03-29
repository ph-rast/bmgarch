# Contributing to `bmgarch`

Contributions to `bmgarch` follow the same principles as stated in the guideline for contributors to the `tidyverse` ecosystem of R (see here for more details: [**development contributing guide**](https://rstd.io/tidy-contrib) ).

## Fixing typos

You can refer to [this](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) for proposing changes with pull request from a fork. 
You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 

If you find any typos in an `.Rd` file, then please make changes in the corresponding `.R` file, as `.Rd` files are automatically generated by [roxygen2](https://roxygen2.r-lib.org/articles/roxygen2.html) and should not be edited by hand.

## Bigger changes

The first step here for you to make any changes is to install `devtools` using `install.packages("devtools")`.
If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal reproducible example, a
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

To contribute a change to `bmgarch`, you follow these steps:

1. Create a branch in git, give it a descriptive name, and make your changes.
2. Push branch to GitHub and issue a pull request (PR).
3. Discuss the pull request.
4. Iterate until either we accept the PR or decide that it's not a good fit for `bmgarch`.

If you're not familiar with git or GitHub, please start by reading
<http://r-pkgs.had.co.nz/git.html>

Pull requests will be evaluated against a checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the need for change. Please describe the problem your PR addresses and show how your pull request solves it as concisely as possible.

Also, include this motivation in `NEWS` so that when a new release of
`bmgarch` comes out it's easy for users to see what has changed. Add your
item at the top of the file and use markdown for formatting. The
news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request, please check to make sure that you haven't accidentally included any unrelated changes. These make it harder to see exactly what's changed, and to evaluate any unexpected side effects.

Each PR corresponds to a git branch, so if you expect to submit multiple changes
make sure to create multiple branches. If you have multiple changes that depend
on each other, start with the first one and don't submit any others until the
first one has been processed.

3. __Use `bmgarch` coding style__. We tried to adhere as closely as possible to the [official `tidyverse` style guide](http://style.tidyverse.org) -- please do so as well. Maintaining a consistent style across the whole code base makes it much easier to jump into the code. If you're modifying existing `bmgarch` code that doesn't follow the style guide, a separate pull request to fix the style would be greatly appreciated.

4. If you're adding new parameters or a new function, you'll also need to document them with [`roxygen2`](https://github.com/klutometis/roxygen). Make sure to re-run `devtools::document()` on the code before submitting.
