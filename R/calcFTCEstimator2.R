calcFTCEstimator2 <-
function(sample, evalBase, evalDer, coefs, small_dom, comp_dom){
	# !!!ATTENTION!!
	# This estimator projects onto the full domain i.e. it extrapolates
	# This estimator is NOT used in the paper NOR in the simulation
	# Exkurs 
	# i <- 10
	# plot(sample[, i], type = "l", lwd = 2)
	# derivative <- evalDer %*% coefs[,i]
	# lines(x = 1:501, integrate_xy(comp_dom,derivative) + sample[1, i], col ="red", lty = 2)
	
	
	# Typical Estimator on usual domain
	mean_sd  			<- rowMeans(sample[1:length(small_dom), ])
	
	#mean_sd2 			<- rowMeans(sample,na.rm = TRUE)
	#kraussMean      	<- apply(sample, 1, mean, na.rm = TRUE)
	#mean_sd[1:10]
	#mean_sd2[1:10]
	#kraussMean[1:10]
	
	
	# FTC Estimator on the rest of the Domain
	X_bar_prime 		<- rowMeans(evalDer %*% coefs)[(length(small_dom)+1):length(comp_dom)]
	
	X_bar_prime_int 	<- integrate_xy(x = comp_dom[(length(small_dom)+1):length(comp_dom)], y = X_bar_prime)
	
	mean_ftc 			<- mean_sd[length(small_dom)] + X_bar_prime_int
	
	# plot
	
	return(c(mean_sd, mean_ftc))
}
