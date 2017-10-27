selBasisLength <-
function(X, comp_dom, base = "fourier", basis_seq = seq(3,51,2)){

	# dim(X) # a p x n matrix
	X_dom   <- apply(X, 2, function(curve) !is.na(curve))
	
	n       <- dim(X)[2]
	# each basis dimension
	bicsPerBasisDim <- lapply(basis_seq, function(nbasis){ # nbasis <- basis_seq[1]
						estBasis 	<- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = nbasis)
						#print(nbasis)
						# each cruve
						bics <- sapply(1:dim(X)[2], function(i){ # i =2
							# print(i)
							argvals <- comp_dom[na.omit(X_dom[, i])]
							y 		<- as.numeric(na.omit(X[, i])) # length(y)
							X_mat 	<- eval.basis(evalarg = comp_dom[na.omit(X_dom[, i])], basisobj = estBasis) #dim(X)
							
							
							# In X there is a constant because of Fourier first basis
							return(bic 	<- BIC(lm(y ~ -1 + X_mat)))
							
							# via FDA Package:
							# estCoef <- smooth.basis(argvals = argvals, y = y, fdParobj = estBasis)
							# estCoef$fd$coefs
							# plot(x=argvals, y = y, type = "l")
							# 
						})
						return(bics)
	})
	names(bicsPerBasisDim) <- basis_seq
	
	# To save computational time, for each b, one has the corresponding bic
	# now on has to rearragne the whole thing such that for each curve, one has all bics
	
	bicMat 	<- matrix(NA, nrow = length(basis_seq), ncol = n)
	rownames(bicMat) <- basis_seq
	for(i in 1:n){
		bicMat[ ,i] <- sapply(bicsPerBasisDim, function(basis) basis[[i]])
	}
	
	return(bicMat)
}
