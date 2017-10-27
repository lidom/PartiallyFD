simDataIC <-
function(n, xi_1_mean, ev, suppliedSettings){ # xi_1_mean = 5
	# Simultes independent continuous data from a normal for x and a uniform for d
	
	settings <- list(sd_part = 1, upper = c(Inf, 1), lower = c(-Inf, .5), xi_1_d_cor = 0.8, d_mean = 0.75, d_var = 1)       
    
    
    matching 				<- intersect(names(settings), names(suppliedSettings))
    settings[matching] 	    <- suppliedSettings[matching]
	
    d 		<- runif(n, min = settings$lower[2], max = settings$upper[2])
	
	x 		<- rnorm( n = n, mean = xi_1_mean, sd = sqrt(ev))
	x 		<- cbind(x,d)
	colnames(x) <- c("x", "d")
	return(sim_data 	= x) # plot(sim_data)
}
