\name{var.class}
\alias{var.class}
\title{Detect the classes of the variables.}
\usage{
var.class(nametable.class, dataset.class)
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
}
\value{
  A vector matching the rows of 'nametable.class'. The
  value includes NA if any variable are only NA's.
}
\description{
  This function gives an initial guess of the classes of
  each variable in the merged data.
}
\examples{
a = data.frame(aa = 1:5, ab = LETTERS[6:2], ac = as.logical(c(0, 1, 0, NA, 0)))
b = data.frame(b1 = letters[12:14], b2 = 3:1)
dat = list(a, b)
name = matrix(c("ab", "aa", "ac", "b1", "b2", NA), ncol = 2)
var.class(name, dat)
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
