\name{plotOrdinalFittedVariable}
\alias{plotOrdinalFittedVariable}
\title{
Function that gives the posibility for the user for plotting in the reduced space an ordered and fitted categorical variable.
}
\description{
Graphical representation of a polytomous ordered variable previously fitted in the reduced space,
according to the Ordinal Logistic Biplot theory. It can be choosen some parameters related to the
way in which the variable is plotted.
}
\usage{
plotOrdinalFittedVariable(nameVariable, coeffic, D,numFactors, planex = 1, planey = 2,
xi = -3.5, xu = 3.5, yi = -3.5, yu = 3.5, margin = 0,
CexVar = 0.7, ColorVar = "blue",
PchVar = 0.7, addToPlot = FALSE, showIIC = TRUE,
iicxi = -2.5, iicxu = 2.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nameVariable}{
Name of the variable the user wants to plot.
}
  \item{coeffic}{
Vector with the estimated coefficients and the thresholds in this order.
}
  \item{D}{
Parameter of the graded response model. In case of coefficients have been estimated by Mirt this parameter should be 1.702.
In other cases it should be 1.
}
  \item{numFactors}{
Number of dimensions of the solution
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
 \item{xi}{
Minimum value on the x-axis.
}
  \item{xu}{
Maximum value on the x-axis.
}
  \item{yi}{
Minimum value on the y-axis.
}
  \item{yu}{
Maximum value on the y-axis.
}
  \item{margin}{
This value stablishes the space between the plotted items and the border of the window.
}
  \item{CexVar}{
Size of the category points. It can be an array with the cex information for each variable.
}
  \item{ColorVar}{
Color for the variables. It can be an array with the color information for each variable.
}
  \item{PchVar}{
Symbol for the variables. It could be an array with the pch information for each variable.
}
  \item{addToPlot}{
Boolean parameter to decide if the user wants to add the ordinal variable
representation to an existing plot.
}
  \item{showIIC}{
Boolean parameter to decide if the user wants to see the item information
curves for each variable. Default value is FALSE.
}
  \item{iicxi}{
Lower limit for the X-axis when plotting item information curves.
}
  \item{iicxu}{
Upper limit for the X-axis when plotting item information curves.
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd,sFormula=NULL,
      numFactors=2,method="EM",penalization=0.2)
nameVariable="Salary"
coeffic = c(olbo$coefs[1,],olbo$thresholds[1,])
plotOrdinalFittedVariable(nameVariable,coeffic,D=1,numFactors = 2)

}

\keyword{plot}