calcFTCEstimator <-
function(fds, evalBase, evalDer, coefs, small_dom, comp_dom, derFds = FALSE){

	# Typical Estimator on usual domain
	mean_sd  			<- rowMeans(fds[1:length(small_dom), ])
	
	# matplot(x = comp_dom, y = fds, type = "l", col = addAlpha("darkgrey", 1), ylim = c(0,5000))
	# lines(x= small_dom, y = mean_sd, col = "darkblue", lwd = 4)
	
	# if the derivatives are not supplied, use the basis which was supplied
	if( isTRUE(!derFds) ) { 
	
		# Take the first derivatives of the sample
		X_prime 			<- evalDer %*% coefs
		
		# Connect the NAs of the sample with ones in the derivatives
		derFds 				<- sapply(1:dim(fds)[2], function(i){ # i = 121
								#print(i)
								# (evalDer %*% coefs)[ , i]
								if(sum(is.na(fds[,i])) != 0){
									v <- c(X_prime[!is.na(fds[,i]), i], rep(NA, times = sum(is.na(fds[,i]))))
									return(v)
								} else {
									return(X_prime[!is.na(fds[,i]), i])
								}
							})
		
		#if(!isTRUE(forceDer)
		
		
		# which.max(abs(derFds[1,]))
		# plot(x = comp_dom, y = derFds[, 183])
		 # dim(coefs)
		 # rev(order(abs(coefs[ ,183])))
			# dates[183]
		
		# FTC Estimator on the rest of the Domain
		X_bar_prime 		<- rowMeans(derFds, na.rm = TRUE)
		
		
	} else { 
		X_bar_prime 		<- rowMeans(derFds, na.rm = TRUE)
	}
	
	# pdf("ableitungen.pdf")
	# matplot(x = comp_dom, y = derFds, type = "l", col = addAlpha("darkgrey", 1), ylim = c(-50,50))
	# lines(x= comp_dom[(length(small_dom)+1):length(comp_dom)], y = X_bar_prime[(length(small_dom)+1):length(comp_dom)], col = "darkblue", lwd = 8)
	# dev.off()
	
	X_bar_prime_int 	<- integrate_xy(x = comp_dom[(length(small_dom)+1):length(comp_dom)], y = X_bar_prime[(length(small_dom)+1):length(comp_dom)])
	
	mean_ftc 			<- mean_sd[length(small_dom)] + X_bar_prime_int
	
	return(list(mean = c(mean_sd, mean_ftc), firstDer = X_bar_prime))
}
