\name{intersect2}
\alias{intersect2}
\title{Obtain the intersection of a list of vectors.
Function "intersect" in the base package can only intersect two
vectors. The function "intersect2" is designed to obtain the
intersection and the difference for more than two vectors. The
input should be a list whose elements are the vectors, and the
outputs include the intersection of all vectors and a list whose
elements are the input vectors substracting the intersection.
Besides, intersect2 allows the labels of the vectors. If a list of
labels is given in the input, then the outputs will also include a
matrix of labels which match the intersection for the vectors, and
a list of labels which match the left part of the vectors.}
\usage{
intersect2(vname, simplifiedname = vname)
}
\arguments{
  \item{vname}{A list of labels.}

  \item{simplifiedname}{A list of vectors to make the
  intersection. Each element in the list has the same
  length as the corresponding element in vname. Default to
  be vname. If simplifiedname is not vname, then it works
  as the real vectors to match, and vname is like the
  labels of simplifiedname. If simplifiedname is the same
  as vname, then the returned value simpleuniq=uniq.}
}
\value{
  The outputs are 'public', 'individual', 'uniq', and
  'simpleuniq'.  'public' is a vector of the intersection
  of 'simplifiedname'.  'individual' is a matrix with the
  original colnames matched to 'public' in all files.
  'simpleuniq' is a list of the left part of
  'simplifiedname' if we pick 'public' out. 'uniq' is a
  list of the left part of 'vname' if we pick 'individual'
  out.
}
\description{
  Obtain the intersection of a list of vectors. Function
  "intersect" in the base package can only intersect two
  vectors. The function "intersect2" is designed to obtain
  the intersection and the difference for more than two
  vectors. The input should be a list whose elements are
  the vectors, and the outputs include the intersection of
  all vectors and a list whose elements are the input
  vectors substracting the intersection. Besides,
  intersect2 allows the labels of the vectors. If a list of
  labels is given in the input, then the outputs will also
  include a matrix of labels which match the intersection
  for the vectors, and a list of labels which match the
  left part of the vectors.
}
\examples{
a = list(x1 = c("label11", "label12"), x2 = c("label21", "label22", "label23"), 
    x3 = c("label31", "label32"))
b = list(x1 = c(1, 2), x2 = c(3, 1, 2), x3 = c(2, 1))
intersect2(a, b)
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
