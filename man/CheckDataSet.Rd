\name{CheckDataSet}
\alias{CheckDataSet}
\title{
Check a data set.
}
\description{
This function checks if a data set is a data frame or a matrix and saves the data as a matrix of integers, and stores the names of rows, colums and levels for each variable as vectors to use them later.
}
\usage{
CheckDataSet(datanom)
}
\arguments{
  \item{datanom}{
It can be a data frame or a matrix.
}
}
\details{
The function checks if some variable has NA values and deletes the corresponding row.
Also checks for missing categories and recodifies the variable keeping the original labels for levels.
}
\value{
An object of class \code{"data.ordinal"}.  This has components:
 \item{datanom }{Matrix of integers with the values of the variables}
 \item{RowNames}{Vector with the names of the rows}
 \item{ColumNames }{Vector with the names of the variables}
 \item{LevelNames }{Levels of each variable}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\examples{
data(LevelSatPhd)
dataChecked = CheckDataSet(LevelSatPhd)
}
\keyword{check}
\keyword{data}
