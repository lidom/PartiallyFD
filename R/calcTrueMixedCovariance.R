calcTrueMixedCovariance <-
function(basis_fd, seq.ev, comp_dom, x_point) { # Lfdobj = 1
	
	# Calculates the mixture covariance rho(X'(s),X(t))
	evalBasis 		<- eval.basis(evalarg = comp_dom, basisobj = basis_fd, Lfdobj=0) # a p x b object
	evalDerBasis 	<- eval.basis(evalarg = comp_dom, basisobj = basis_fd, Lfdobj=1) # a p x b object
	
	singleCovMats <- lapply(1:length(seq.ev), function(ind){ # ind = 1
								seq.ev[ind] * (evalBasis[x_point ,ind] %*% t(evalDerBasis[ ,ind]))
							})
	# Calculate true covariance via cov(s,t) = sum_i=1^b lambda_i * psi_i(t)*psi_i(s)
	trueCov   <- singleCovMats[[1]]
	
	for(i in 2:length(seq.ev)){ # i = 2; isSymmetric(singleCovMats[[i]])
		trueCov 	<- trueCov + singleCovMats[[i]]
	}
	
	# library(rgl)
	# persp3d(x = comp_dom, y =comp_dom, z = trueCov, col="skyblue")# isSymmetric(cov_hat_III)
	
	return(trueMixCov = trueCov)
	
}
