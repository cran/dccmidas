\name{Det}
\alias{Det}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix determinant
}
\description{
Calculates the determinant of a numeric matrix.
}
\usage{
Det(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric matrix
}
}
\value{
The determinant of x.
}

\examples{
x<-matrix(sample(1:25,25,replace=TRUE),ncol=5)
Det(x)
}