\name{scale_kstest}
\alias{scale_kstest}
\title{Compute the p-values of the Kolmogorov-Smirnov tests between
different sources for each variable.
This function is used to detect whether the matched variables from
different files have different distributions. For each variable,
it will compute the pairwise KS-test p-values among the sources,
then report the lowest p-value as the indice for this variable.}
\usage{
scale_kstest(nametable.class, dataset.class, name.class, varclass = NULL)
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

  \item{varclass}{A character vector of variable classes.
  The length of the vector must be equal to the number of
  rows in nametable.class. All the classes should be in
  "numeric", "integer", "factor", and "character". Default
  to be null, then it will be determined by
  \code{\link{var.class}}.}
}
\value{
  A vector of p-values from the KS-test for each
  variable.The p-values are between 0 and 1, or equal to 9
  if one of more groups only have NA's.
}
\description{
  Compute the p-values of the Kolmogorov-Smirnov tests
  between different sources for each variable. This
  function is used to detect whether the matched variables
  from different files have different distributions. For
  each variable, it will compute the pairwise KS-test
  p-values among the sources, then report the lowest
  p-value as the indice for this variable.
}
\examples{
a = data.frame(aa = 1:5, ab = LETTERS[6:2], ac = as.logical(c(0, 1, 0, NA, 0)))
b = data.frame(b1 = letters[12:14], b2 = 3:1)
dat = list(a, b)
name = matrix(c("ab", "aa", "ac", "b1", "b2", NA), ncol = 2)
colnames(name) = c("a", "b")
newname = c("letter", "int", "logic")
scale_kstest(name, dat, newname)
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
