\name{BiplotDensity}
\alias{BiplotDensity}
\title{
Density plot of a variable with contours overlaid.
}
\description{
This function draws for a set of points a density contour lines biplot. This tool uses the kde2d function of library MASS.
}
\usage{
BiplotDensity(X, y = NULL, nlevels = max(y), grouplabels = 1:nlevels,
  ncontours = 6, groupcols = 1:nlevels, img = TRUE, separate = FALSE,
  ncolors = 20, ColorType = 4, xliml = -1, xlimu = 1, yliml = -1,
  ylimu = 1, plotInd = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix with the items coordinates for the plane in which the the contour lines will be plotted.
}
  \item{y}{
Categorical variable that the user wants to study.
}
  \item{nlevels}{
Maximum value of the variable specified in the second parameter.
}
  \item{grouplabels}{
Set of labels for the centers of each category. It should be a vector.
}
  \item{ncontours}{
Number of contours that will be used in the representation.
}
  \item{groupcols}{
Set of colors for each contour lines group. It should be a vector.
}
  \item{img}{
Boolean parameter to specify if an image will be plotted behind the contour lines indicating with diferent colorations the density of the points.Default value is TRUE. 
}
  \item{separate}{
Boolean parameter to specify for each category of the variable a separate biplot representation.Default value is FALSE. 
}
  \item{ncolors}{
Number of colors that will be used in the image.
}
  \item{ColorType}{
Type of coloration for the image in the biplot. It should be a number between 1 and 5.
}
 \item{xliml}{
Minimum value on the x-axis.
}
  \item{xlimu}{
Maximum value on the x-axis.
}
  \item{yliml}{
Minimum value on the y-axis.
}
  \item{ylimu}{
Maximum value on the y-axis.
}
  \item{plotInd}{
Should the individuals be plotted?
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd[,1:4])
x = olbo$ItemsCoords[, 1]
y = olbo$ItemsCoords[, 2]
plot(x,y, cex = 0, xlim=c(-1,1),ylim=c(-1,1))
X = olbo$ItemsCoords
y = as.matrix(as.numeric(LevelSatPhd[,4]))
gcols = c("midnightblue","black","red","gray87")
BiplotDensity(X,y,groupcols = gcols)
}
\keyword{Density}