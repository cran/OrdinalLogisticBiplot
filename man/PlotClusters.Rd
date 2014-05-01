\name{PlotClusters}
\alias{PlotClusters}
\title{
Graphical representation of clusters of individuals.
}
\description{
This function uses a nominal variable to represent groups or clusters of individuals. The clusters can be the result of a clustering algorithm or the groups defined by a external nominal variable. The centroids and convex hulls for each cluster can be represented.

}
\usage{
PlotClusters(A, Groups = ones(c(nrow(A), 1)),
        colors = NULL, chulls = TRUE, centers = TRUE, ConfidentLevel = 0.95)
}
\arguments{
  \item{A}{
A matrix with the coordinates of each point. It should have only two columns.
}
  \item{Groups}{
Clustering variable: the cluster for each observation.
}
  \item{colors}{
It is a vector used to specify the color for each cluster. 
}
  \item{chulls}{
Should convex hulls regions for each cluster be plotted?
}
  \item{centers}{
Should centroids of each cluster be plotted?
}
  \item{ConfidentLevel}{
Numerical value between 0 and 1. If it's value is 0.95, five percent of the points with higher distances
to the center of each cluster will not be used to calculate centroids and convex hulls.
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd)
x = olbo$RowCoords[, 1]
y = olbo$RowCoords[, 2]
plot(x,y, cex = 0.8, pch=17, xlim=c(-2,2),ylim=c(-2,2))
GroupsF = as.factor(LevelSatPhd[,4])
PlotClusters(olbo$RowCoords, Groups = GroupsF,
    colors = c(1,2,3,4),chulls = TRUE,centers = TRUE,ConfidentLevel=NULL)

}
\keyword{cluster}
