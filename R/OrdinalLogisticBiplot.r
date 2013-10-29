# file OrdinalLogisticBiplot/R/OrdinalLogisticBiplot.R
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


OrdinalLogisticBiplot <- function(datanom,sFormula=NULL,numFactors=2,coordinates="EM",rotation="varimax",
                          metfsco="EAP",nnodos = 10, tol = 1e-04, maxiter = 100, penalization = 0.1,cte=TRUE,
                          show=FALSE,ItemCurves = FALSE){

  #We have to check if datanom is a data frame or a matrix to tackle with sFormula parameter
  if(!(is.null(sFormula))){
    if(is.data.frame(datanom)){
      datanom = model.frame(formula=sFormula,data=datanom)
    }else if(is.matrix(datanom)){
              datanom = model.frame(formula=sFormula,data=as.data.frame(datanom))
              datanom = as.matrix(datanom)
          }else{
            print("It is not posible to use the formula passed as parameter. Data are not a data frame nor matrix")
          }
  }
  if(ncol(datanom) <= numFactors){
    stop("It is not posible to reduce dimension because number of Factors is lower than variables in our data set")
  }
  dataSet = CheckDataSet(datanom)
  datanom = dataSet$datanom

  nRowsdata <- dim(datanom)[1]
  nColsdata <- dim(datanom)[2]
  numVar <- ncol(datanom)    
  datanomcats = apply(datanom[,1:numVar], 2, function(x) nlevels(as.factor(x)))

  x <- matrix(0,nRowsdata,numFactors)
  
  if (coordinates == "EM"){
    olb = OrdinalLogBiplotEM(datanom,dim = numFactors, nnodos = nnodos, tol = tol, maxiter = maxiter, penalization = penalization)
    catOrdBiplotPenal = GetOrdinalBiplotObjectPenal(dataSet$ColumNames,olb)
    coefs = catOrdBiplotPenal$models$ColumnParameters$coefficients
    thresholds = catOrdBiplotPenal$models$ColumnParameters$thresholds
    x = catOrdBiplotPenal$models$RowCoordinates
    VariableModels = catOrdBiplotPenal$matBiplot
    Fitting = catOrdBiplotPenal$models$ColumnParameters$fit
    LogLik = catOrdBiplotPenal$models$LogLikelihood
    rotation="Not applicable"
    metfsco="Not applicable"

   }else if (coordinates == "MIRT"){           
      rows = EstimationRowsMIRT(datanom,numFactors = numFactors,metfsco=metfsco,rotation=rotation)  #Se estima una vez y luego se elige el plano que queremos estudiar en GetOrdinalBiplotObjectMIRT
      catOrdBiplot = GetOrdinalBiplotObjectMIRT(rows)        
      x=catOrdBiplot$scores                   
      VariableModels = catOrdBiplot$matBiplot
      coefs = catOrdBiplot$sepCoefMirt$xCoeffic
      thresholds = catOrdBiplot$sepCoefMirt$indCoeffic
      Fitting = NULL
      LogLik = NULL
      rotation=rotation
      metfsco=metfsco
      nnodos="Not applicable"
  }
  
   dimnames(x)[[1]]=dimnames(datanom)[[1]]    #c(1:nRowsdata)
   dimnames(x)[[2]]=c(1:numFactors)

    ordinal.logistic.biplot<-list()
    ordinal.logistic.biplot$dataSet = dataSet
    ordinal.logistic.biplot$ItemsCoords = x
    ordinal.logistic.biplot$VariableModels = VariableModels
    ordinal.logistic.biplot$Fitting = Fitting 
    ordinal.logistic.biplot$coefs = coefs  
    ordinal.logistic.biplot$thresholds = thresholds   
    ordinal.logistic.biplot$NumFactors = numFactors
    ordinal.logistic.biplot$Coordinates = coordinates
    ordinal.logistic.biplot$Rotation = rotation
    ordinal.logistic.biplot$Methodfscores = metfsco
    ordinal.logistic.biplot$NumNodos = nnodos
    ordinal.logistic.biplot$tol = tol
    ordinal.logistic.biplot$maxiter = maxiter
    ordinal.logistic.biplot$penalization = penalization
    ordinal.logistic.biplot$cte = cte
    ordinal.logistic.biplot$show = show
    ordinal.logistic.biplot$ItemCurves =  ItemCurves
    ordinal.logistic.biplot$LogLik =  LogLik

    class(ordinal.logistic.biplot)='ordinal.logistic.biplot'
    return(ordinal.logistic.biplot)
}

summary.ordinal.logistic.biplot <- function(object,...) {
      x = object

      if(x$Coordinates == "EM"){
         	cat(paste("Ordinal Logistic Biplot Estimation", "with Ridge Penalizations", x$penalization, " and logit link"), "\n")        
         	cat("logLikelihood: ", x$LogLik, "\n")
         	cat("\n\nPercentage of correct classifications and Pseudo R-squared measures: \n")
         	print(x$Fitting)
      }else{
         	cat(paste("Ordinal Logistic Biplot Estimation", "with MIRT method"), "\n")
          cat(paste("Rotation", x$Rotation), "\n")
          cat(paste("Method of fscores calculation:", x$Methodfscores), "\n")    
      }
    	cat("n: ", nrow(x$dataSet$datanom), "\n")
    	cat("Coordinates for the individuals: ","\n")
      print(x$ItemsCoords)
     	cat("\n Coefficients:\n")
     	print(x$coefs)
     	cat("\n Thresholds:\n")
     	print(x$thresholds)
}

plot.ordinal.logistic.biplot <- function(x,planex=1,planey=2,AtLeastR2 = 0.01,
        xlimi=-1.5,xlimu=1.5,ylimi=-1.5,ylimu=1.5,margin = 0,
        ShowAxis = TRUE, PlotVars = TRUE, PlotInd = TRUE, LabelVar = TRUE,
        LabelInd = TRUE,CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL,
        PchInd = NULL, PchVar = NULL,showIIC=FALSE,iicxi=-1.5,iicxu=1.5,legendPlot = FALSE,...) {
  olbo = x
  n = nrow(olbo$dataSet$datanom)
	p = ncol(olbo$dataSet$datanom)
	RowNames = olbo$dataSet$RowNames
	VarNames = olbo$dataSet$ColumNames

	DimNames = "Dim 1"
	for (i in 2:olbo$NumFactors)
     DimNames = c(DimNames, paste("Dim", i))

  # Determining sizes and colors of the points
	if (is.null(CexInd)){
		CexInd = rep(0.5, n)
	}else{
     if (length(CexInd) == 1){
       CexInd = rep(CexInd, n)
     }else if(length(CexInd) < n){
             CexInd = rep(0.5, n)
           }else{
             CexInd = CexInd[1:n]
           }
	}

  if (is.null(ColorInd)){
  		ColorInd = rep("black", n)
	}else{
     if (length(ColorInd) == 1){
       ColorInd = rep(ColorInd, n)
     }else if(length(ColorInd) < n){
             ColorInd = rep("black", n)
           }else{
             ColorInd = ColorInd[1:n]
           }
	}

 	if (is.null(PchInd)){
		PchInd = rep(1, n)
	}else{
     if (length(PchInd) == 1){
       PchInd = rep(PchInd, n)
     }else if(length(PchInd) < n){
             PchInd = rep(1, n)
           }else{
             PchInd = PchInd[1:n]
           }
	}

	if (is.null(CexVar)){
		CexVar = rep(0.8, p)
	}else{
     if (length(CexVar) == 1){
       CexVar = rep(CexVar, p)
     }else if(length(CexVar) < p){
             print("It has been specified lower cex values for the variables than variables")
             CexVar = rep(0.8, p)
           }else{
             CexVar = CexVar[1:p]
           }
  }

	if (is.null(PchVar)){
    PchVar = c(0:(p-1))
	}else{
    if (length(PchVar) == 1){
       PchVar = c(0:(p-1))
     }else if(length(PchVar) < p){
             print("It has been specified lower pch values for the variables than variables")
             PchVar = c(0:(p-1))
           }else{
             PchVar = PchVar[1:p]
           }
  }
  if (is.null(ColorVar)){
		ColorVar = colors()[20 + 2*c(1:p)]
	}else{
     if (length(ColorVar) == 1){
       ColorVar = colors()[20 + 2*c(1:p)]
     }else if(length(ColorVar) < p){
             print("It has been specified lower color values for the variables than variables")
             ColorVar = colors()[20 + 2*c(1:p)]
           }else{
             ColorVar = ColorVar[1:p]
           }
  }

	if (ShowAxis) {
		xaxt = "s"
		yaxt = "s"
	} else {
		xaxt = "n"
		yaxt = "n"
	}

  if ((margin < 0) | (margin > 0.3))
		margin = 0

  xmin= xlimi - (xlimu - xlimi) * margin
  xmax= xlimu + (xlimu - xlimi) * margin
  ymin= ylimi - (ylimu - ylimi) * margin
  ymax= ylimu + (ylimu - ylimi) * margin

  dev.new()
  if(PlotInd == TRUE){
    plot(olbo$ItemsCoords[,planex], olbo$ItemsCoords[,planey], cex = CexInd, col=ColorInd, pch = PchInd, asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        main="Ordinal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))
    if(LabelInd == TRUE){
         text(olbo$ItemsCoords[,planex], olbo$ItemsCoords[,planey],row.names(RowNames), cex = CexInd,col=ColorInd,pos=1,offset=0.1)
    }
  }else{
     plot(olbo$ItemsCoords[,planex], olbo$ItemsCoords[,planey], cex = 0,asp=1, xaxt = xaxt, yaxt = yaxt ,xlim=c(xmin,xmax),ylim=c(ymin,ymax),
        main="Ordinal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))
  }

  if(PlotVars){
      if(olbo$Coordinates == "EM"){
        D = 1
      }else{
        D = 1.702
      }
      plot.ordinalBiplotPenal(olbo,D,planex,planey,xi=xlimi,xu=xlimu,
              yi=ylimi,yu=ylimu,margin = margin, CexVar = CexVar,ColorVar = ColorVar,
              PchVar = PchVar)
      if(legendPlot){
        legend("bottomright", legend=VarNames, col= ColorVar,pch=PchVar,cex=0.7)
      }
  }
  
  if(showIIC){
    for(v in 1:ncol(olbo$dataSet$datanom)){
      nameVariable = VarNames[v]
      coeffic = olbo$VariableModels[,v]$coef
      slopeort = olbo$VariableModels[,v]$slope
      D = 1
      numcat = max(olbo$dataSet$datanom[,v])
      xi = iicxi
      xu = iicxu
      if(olbo$Coordinates == "MIRT"){
        D = 1.702
      }
      plotCurvesCategoriesVariable(coeffic,slopeort,D,numcat,nameVariable,xi,xu,planex,planey)
    }
  }
  
}
