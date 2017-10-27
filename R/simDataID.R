simDataID <-
function(n, xi_1_mean, ev, suppliedSettings){ # xi_1_mean = 5
	# Simulates independent discrete data from a n distribution
	
	# default settings
    settings <- list(sd_part = 1, upper = c(Inf, 1), lower = c(-Inf, .5), xi_1_d_cor = 0.8, d_mean = 0.75, d_var = 1)       
    

    matching 				<- intersect(names(settings), names(suppliedSettings))
    settings[matching] 	    <- suppliedSettings[matching]

	simData_id 			<- cbind(	x 	<- rnorm(n = n, mean = xi_1_mean, sd = sqrt(ev)),
									d	<- sample( x = c(0.5, 1), 		size = n, replace = TRUE))
	colnames(simData_id) 	<- c("x", "d")
	return(sim_data 	= simData_id)
}
