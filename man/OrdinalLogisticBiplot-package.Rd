\name{OrdinalLogisticBiplot-package}
\alias{OrdinalLogisticBiplot-package}
\docType{package}
\title{
Ordinal Logistic Biplot representations for polytomous ordered data.
}
\description{
Analysis of a matrix of polytomous ordered items using Ordinal Logistic Biplots (OLB). The OLB procedure extends the binary logistic biplot to ordinal (polytomous) data. 
}
\details{
\tabular{ll}{
Package: \tab OrdinalLogisticBiplot\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-10-29\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon
Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\keyword{package}
\seealso{
\code{\link{OrdinalLogisticBiplot}},\code{\link{OrdinalLogBiplotEM}}
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd,sFormula=NULL,numFactors=2,
coordinates="EM",penalization=0.2,show=FALSE)
summary(olbo)
plot(olbo,PlotInd=TRUE,xlimi=-1,xlimu=1,ylimi=-1,ylimu=1,margin = 0.2,
ColorVar = c("red","green","black","blue","yellow"),CexVar = c(0.7),showIIC=FALSE)

}
