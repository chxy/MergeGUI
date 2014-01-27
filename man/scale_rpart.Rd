\name{scale_rpart}
\alias{scale_rpart}
\title{Compute the misclassification rate for each variable.
When merging data from several datasets, it is meaningful to
detect whether the matched variables from different files have
different centers. The function computes the misclassification
rate variable by variable using classification tree (the rpart
package). It will firstly merge the dataset by the given
nametable.class, then use rpart for each variable to seperate the
data without any covariates and compute the misclassification
rate.}
\usage{
scale_rpart(nametable.class, dataset.class, name.class, varclass = NULL)
}
\arguments{
  \item{nametable.class}{A matrix of the matched variable
  names. The number of columns is equal to the number of
  files. The column names are required.  Each row
  represents a variable that is going to be merged. Any
  elements except NA in nametable.class must be the
  variable names in dataset.class.}

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
  A vector of the misclassification rate. The rate is
  between 0 and 1, or equal to 9 if one of more groups only
  have NA's.
}
\description{
  Compute the misclassification rate for each variable.
  When merging data from several datasets, it is meaningful
  to detect whether the matched variables from different
  files have different centers. The function computes the
  misclassification rate variable by variable using
  classification tree (the rpart package). It will firstly
  merge the dataset by the given nametable.class, then use
  rpart for each variable to seperate the data without any
  covariates and compute the misclassification rate.
}
\examples{
a = data.frame(aa = 1:5, ab = LETTERS[6:2], ac = as.logical(c(0, 1, 0, NA, 0)))
b = data.frame(b1 = letters[12:14], b2 = 3:1)
dat = list(a, b)
name = matrix(c("ab", "aa", "ac", "b1", "b2", NA), ncol = 2)
colnames(name) = c("a", "b")
newname = c("letter", "int", "logic")
scale_rpart(name, dat, newname)
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
