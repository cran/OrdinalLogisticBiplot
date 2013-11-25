# file OrdinalLogisticBiplot/R/EstimationRowsMIRT.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

EstimationRowsMIRT <- function(dataFile,numFactors=2,metfsco="EAP",rotation="varimax",maxiter=100)
{
  nRowsdataF <- dim(dataFile)[1]
  nColsdataF <- dim(dataFile)[2]
  varStudy = matrix(0,nRowsdataF,nColsdataF)
  for(i in 1:nColsdataF){
    varStudy[,i]= factor(dataFile[,i])
  }
  dimnames(varStudy)[[2]]=dimnames(dataFile)[[2]]
  datanomcats = apply(varStudy, 2, function(x) nlevels(as.factor(x)))

  technical = list()
  technical$NCYCLES = maxiter

  (modF<-mirt(varStudy,numFactors,rotate=rotation,,technical=technical))
  summ = summary(modF)
  tabscores<-fscores(modF,full.scores=TRUE,method=metfsco)
  xSubi = tabscores[,(nColsdataF+1):(nColsdataF+numFactors)]
  catsVarMax = max(datanomcats)
  sepCoefMirt = ExtractMirtCoefficients(catsVarMax,cmodF = coef(modF),numFactors = numFactors)

  result = list()
  result$estimRows = xSubi
	result$dataFactor = varStudy
	result$rotation = rotation
	result$metfsco = metfsco
	result$numFactors = numFactors
	result$coefMirt = coef(modF)
	result$sepCoefMirt = sepCoefMirt
  result$summ = summ
	class(result) <- "EstimationRowsMIRT"
	return(result)

}
