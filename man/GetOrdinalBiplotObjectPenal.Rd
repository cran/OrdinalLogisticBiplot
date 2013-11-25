\name{GetOrdinalBiplotObjectPenal}
\alias{GetOrdinalBiplotObjectPenal}
\title{
Calculation of the ordinal biplot object.
}
\description{
This function uses the estimation from the EM alternated algorithm to calculate all the information needed for representing the ordinal biplot, as the biplot axis for each variable and the points separating the visible categories.
}
\usage{
GetOrdinalBiplotObjectPenal(ColumNames, olb, planex = 1, planey = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ColumNames}{
Vector with the names of the studied variables.
}
  \item{olb}{
Object of class "ordinal.logistic.biplotEM" with the alternated algorithm estimation.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
}
\details{
This function should only be used if the user selects "EM" coordinates as the estimation method in OrdinalLogisticBiplot function.
}
\value{
An object of class \code{"CategOrdBiplotPen"}.  This has components:
 \item{models}{Estimated object with the EM alternated algoritm of the class "ordinal.logistic.biplotEM". This object has the following items: RowCoordinates(Coordinates for the rows in the reduced space),ColumnParameters(estimated values for the coefficients and the thresholds and indicators of the fitting, like logLik, Deviance, p-values and pseudo R-squared values),loadings(factor loadings values),LogLikelihood(total logarithm of the likelihood),r2(communalities) and Ncats(number of categories of each variable)}
 \item{matBiplot }{Matrix with so many columns as variables. Each column keeps a list of variables describing the biplot elements.The variables are: var(name of the variable), cosines(cosines of the angles between the straight line and the axes), numcat(Number of categories of the variable), coef(estimated coefficients for this variable), slope(slope of the biplot axis), order(boolean variable that is TRUE if the coordinates of the curves intersections [1-2, 2-3, ..., (numcat-1)-numcat] with the biplot axis are ordered), pointsc(matrix with the intersection points between curves calculated as before), pointprob(matrix with the intersection points of each curve with the rest of them and with the biplot axis. It has 5 columns: the first one is the probability for the points in witch the curves are equal, the second an third are the coordinates x and y of the intersection of the two curves with the biplot axis, and the last two are the categories compared. This matrix has zeros if the variable order is TRUE, because it has not been calculated because it is not necessary.)}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{OrdinalLogBiplotEM}},\code{\link{OrdinalLogisticBiplot}}
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd)
ColumNames = olbo$dataSet$ColumNames
olb = olbo$estimObject
catOrdBiplotPenal = GetOrdinalBiplotObjectPenal(ColumNames,olb)
}
\keyword{BiplotObject}
