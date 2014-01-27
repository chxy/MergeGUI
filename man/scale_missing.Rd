\name{scale_missing}
\alias{scale_missing}
\title{Chi-square tests for the counts of missing and non-missing.
This function is used to detect whether the matched variables from
different files have different missing patterns. For each
variable, it will firstly count the missing and non-missing values
among the sources, and then form a contingency table. The p-value
of Chi-square test is computed from the contingency table and
finally reported for the variable.}
\usage{
scale_missing(nametable.class, dataset.class, name.class)
}
\arguments{
  \item{nametable.class}{A matrix of the matched variable
  names. The number of columns is equal to the number of
  files. Each row represents a variable that is going to be
  merged. Any elements except NA in nametable.class must be
  the variable names in dataset.class.}

  \item{dataset.class}{The dataset list. The length of the
  list is equal to the number of files, and the order of
  the list is the same as the order of columns in
  nametable.class.}

  \item{name.class}{A character vector of variable names.
  The length of the vector must be equal to the number of
  rows in nametable.class. Since the variable names in
  nametable.class may not be consistent, name.class is
  needed to name the variables.}
}
\value{
  A vector of p-values from the Chisquare-test for the
  missings of each variable. The p-values are between 0 and
  1.
}
\description{
  Chi-square tests for the counts of missing and
  non-missing. This function is used to detect whether the
  matched variables from different files have different
  missing patterns. For each variable, it will firstly
  count the missing and non-missing values among the
  sources, and then form a contingency table. The p-value
  of Chi-square test is computed from the contingency table
  and finally reported for the variable.
}
\examples{
a = data.frame(aa = 1:5, ab = LETTERS[6:2], ac = as.logical(c(0, 1, 0, NA, 0)))
b = data.frame(b1 = letters[12:14], b2 = 3:1)
dat = list(a, b)
name = matrix(c("ab", "aa", "ac", "b1", "b2", NA), ncol = 2)
colnames(name) = c("a", "b")
newname = c("letter", "int", "logic")
scale_missing(name, dat, newname)
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
