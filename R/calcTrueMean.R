calcTrueMean <- function(basis_fd, xi_means, comp_dom) {
	# Calculates the true mean function for given basis and the statistical properties of the coefficients
	evalBasis <- eval.basis(evalarg = comp_dom, basisobj = basis_fd) # a p x b object
	
	singleMeansMat <- lapply(1:length(xi_means), function(ind){ # ind = 1
								return(xi_means[ind] * evalBasis[ ,ind])
							})
	# Calculate true covariance via cov(s,t) = sum_i=1^b lambda_i * psi_i(t)*psi_i(s)
	trueMean   <- singleMeansMat[[1]]
	
	for(i in 2:length(xi_means)){ # i = 2
		trueMean 	<- trueMean + singleMeansMat[[i]]
	}
	
	# library(rgl)
	# curve(singleMeansMat[[2]]
	# persp3d(x = comp_dom, y =comp_dom, z = trueCov, col="skyblue")# isSymmetric(cov_hat_III)
	
	return(trueMean = trueMean)
	
}
