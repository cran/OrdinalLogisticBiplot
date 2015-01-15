\name{summary.pordlogist}
\alias{summary.pordlogist}
\title{
Summary Method Function for Objects of Class 'pordlogist'
}
\description{
This function shows a summary of the principal results for the estimation
  for individuals and variables, like number of observations, the number of 
  iterations, the covariances matrix, some Pseudo R-squared indices with the
  correct classification percentage of each regression and the logLikelihood with
  "Estimate coefficients", "Std. Error", "z value" or "Pr(>|z|)" values.
}
\usage{
\method{summary}{pordlogist}(object,\dots)
}
\arguments{
  \item{object}{
This parameter keeps 'pordlogist' object for a variable.
}
  \item{\dots}{
Additional parameters to summary.
}
}
\details{
This function is a method for the generic function summary() for class "pordlogist".
It can be invoked by calling summary(x) for an object x of the appropriate class.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{pordlogist}},\code{\link{CheckDataSet}},\code{\link{OrdinalLogBiplotEM}}
}
\examples{
data(LevelSatPhd)
dataSet = CheckDataSet(LevelSatPhd)
datanom = dataSet$datanom
olb = OrdinalLogBiplotEM(datanom,dim = 2, nnodos = 10, tol = 0.001,
  maxiter = 100, penalization = 0.2)
model = pordlogist(datanom[, 1], olb$RowCoordinates, tol = 0.001,
  maxiter = 100, penalization = 0.2)
summary(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary}

