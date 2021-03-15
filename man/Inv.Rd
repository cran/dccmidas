\name{Inv}
\alias{Inv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse of a matrix
}
\description{
Calculates the inverse of a numeric matrix
}
\usage{
Inv(x)
}
\arguments{
  \item{x}{
a numeric matrix
}
}
\value{
The inverse of x.
}

\examples{
x<-matrix(sample(1:25,25,replace=TRUE),ncol=5)
Inv(x)
}