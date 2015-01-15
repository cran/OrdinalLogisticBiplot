\name{plotOrdinalVariable}
\alias{plotOrdinalVariable}
\title{
This function plots in the reduced space an ordered categorical variable.
}
\description{
Graphical representation of a polytomous ordered variable in the reduced space,
according to Vicente-Villardon & Hernandez-Sanchez(2014) methodology.
It can be choosen some parameters related to the way in which the variable is plotted.
}
\usage{
plotOrdinalVariable(ordinalfVar,nameVariable,estimRows,planex = 1,planey = 2,
  xi=-3.5,xu=3.5,yi=-3.5,yu=3.5,margin=0,CexVar=0.7,ColorVar="blue",
  PchVar=0.7,addToPlot=FALSE,showIIC = TRUE,iicxi=-2.5,iicxu=2.5,
  tol = 1e-04, maxiter = 100, penalization = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordinalfVar}{
The ordinal variable. It must be an ordered factor.
}
  \item{nameVariable}{
Name of the variable that the user wants to represent.
}
  \item{estimRows}{
Matrix with the estimated coordinates for the individuals in the reduced dimension.
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
Boolean parameter to decide if the user wants to add the ordinal variable representation to an existing plot.
}
  \item{showIIC}{
Boolean parameter to decide if the user wants to see the item information curves for each variable. Default value is FALSE.
}
  \item{iicxi}{
Lower limit for the X-axis when plotting item information curves.
}
  \item{iicxu}{
Upper limit for the X-axis when plotting item information curves.
}
  \item{tol}{
Tolerance for the iterations.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{penalization}{
Penalization used to avoid singularities.
}
}
\references{
Vicente-Villardon, J. L. & Hernandez, J. C. & (2014) Logistic Biplots for ordinal
data with an application to job satisfaction of doctorate degree holders in Spain.
Preprint available at arXiv.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd,sFormula=NULL,
  numFactors=2,method="EM")
ordinalfVar = factor(LevelSatPhd[,1],ordered=TRUE)
levels(ordinalfVar) = c("VS","SS","SD","VD")
estimRows = olbo$RowCoords
nameVariable = "Salary"   
plotOrdinalVariable(ordinalfVar,nameVariable,estimRows,planex = 1,
  planey = 2,xi=-1.5,xu=1.5,yi=-1.5,yu=1.5,
  margin=0.2,CexVar=0.7,showIIC = TRUE)     
}

\keyword{plot}
