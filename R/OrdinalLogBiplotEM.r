# file OrdinalLogisticBiplot/R/OrdinalLogBiplotEM.R
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


OrdinalLogBiplotEM <- function(x,dim = 2, nnodos = 15, tol = 0.001, maxiter = 100, penalization = 0.2) {
	Q = multiquad(nnodos, dim)
	X = Q$X
	A = Q$A
	q = dim(X)[1]
	p = dim(x)[2]
	Ncats = ColMax(x)
	Maxcat = max(Ncats)
	G = NominalMatrix2Binary(x)
	s = dim(G)[1]
	n = dim(G)[2]
	par = list()
	par$coefficients = array(0, c(p, dim))
	par$thresholds = array(0, c(p, Maxcat - 1))
  par$fit = array(0,c(p,4))                 
  dimnames(par$fit)[[1]]= dimnames(x)[[2]][1:dim(x)[2]]       
  dimnames(par$fit)[[2]]= c("PCC","CoxSnell","Macfaden","Nagelkerke")  
	corr = afc(G, neje = dim)

	ability = corr$RowCoordinates[, 1:dim]

	logLikold = 0
	for (j in 1:p) {
		model = pordlogist(x[, j], ability, tol = tol, maxiter = maxiter, penalization = penalization)
		par$coefficients[j, ] = model$coefficients
    par$thresholds[j,1:nrow(model$thresholds)] = model$thresholds  
		par$fit[j,1] = model$PercentClasif      
		par$fit[j,2] = model$CoxSnell           
		par$fit[j,3] = model$MacFaden         
		par$fit[j,4] = model$Nagelkerke           
		logLikold = logLikold + model$logLik
	}
	error = 1
	iter = 0
	while ((error > tol) & (iter < maxiter)) {
		# E-step - ability estimation
		iter = iter + 1
		PT = EvalOrdlogist(X, par, Ncats)
		L = matrix(1, s, q)
		for (l in 1:s) for (k in 1:q) L[l, k] = prod(PT[k, ]^G[l, ])
		Pl = L %*% A
		ability = matrix(0, s, dim)
		for (l in 1:s) for (j in 1:dim) {
			for (k in 1:q) ability[l, j] = ability[l, j] + X[k, j] * L[l, k] * A[k]
			ability[l, j] = ability[l, j]/Pl[l]
		}
		# M-step  -  Parameter estimation
		logLik = 0
		for (j in 1:p) {
			model = pordlogist(x[, j], ability, tol = tol, maxiter = maxiter, penalization = penalization)
			par$coefficients[j, ] = model$coefficients
      par$thresholds[j,1:nrow(model$thresholds)] = model$thresholds  
  		par$fit[j,1] = model$PercentClasif      
  		par$fit[j,2] = model$CoxSnell           
  		par$fit[j,3] = model$MacFaden         
  		par$fit[j,4] = model$Nagelkerke           
			logLik = logLik + model$logLik
		}
		error = abs((logLik - logLikold)/logLik)
		logLikold = logLik
	}
	d = sqrt(rowSums(par$coefficients^2) + 1)
	loadings = solve(diag(d)) %*% par$coefficients
	thresholds = solve(diag(d)) %*% par$thresholds
	r2 = rowSums(loadings^2)
	model = list()
	model$RowCoordinates = ability
	model$ColumnParameters = par
	model$loadings = loadings
	model$LogLikelihood = logLik
	model$r2 = r2
	model$Ncats=Ncats
	class(model) = "ordinal.logistic.biplotEM"
	return(model)
}

