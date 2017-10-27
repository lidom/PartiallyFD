calcTrueCovariance <-
function(basis_fd, seq.ev, comp_dom, Lfdobj=0) { # Lfdobj = 1
	# Calculates the true covariance function for a basis with eigenvalues on a complete domain.
	# Can calculate the covariance function for the derivatves at this basis as well
	evalBasis <- eval.basis(evalarg = comp_dom, basisobj = basis_fd, Lfdobj=Lfdobj) # a p x b object
	
	singleCovMats <- lapply(1:length(seq.ev), function(ind){ # ind = 1
								seq.ev[ind] * (evalBasis[ ,ind] %*% t(evalBasis[ ,ind]))
							})
	# Calculate true covariance via cov(s,t) = sum_i=1^b lambda_i * psi_i(t)*psi_i(s)
	trueCov   <- singleCovMats[[1]]
	
	for(i in 2:length(seq.ev)){ isSymmetric(singleCovMats[[i]])
		trueCov 	<- trueCov + singleCovMats[[i]]
	}
	
	# library(rgl)
	# persp3d(x = comp_dom, y =comp_dom, z = trueCov, col="skyblue")# isSymmetric(cov_hat_III)
	
	return(trueCov = trueCov)
	
}
