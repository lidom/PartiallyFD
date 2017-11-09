simDataSC <-
function(n, xi_1_mean, ev, suppliedSettings){ # xi_1_mean = 5; n = 1000
# Simulates dependency as standard transformation:
# Let X ~ f_X, then F(X) ~ UNIF[0,1]
# Where UNIF[0,1] will be transformed to UNIF(0.5, 1) and half of them to be 0.5
    
    # default settings
    settings <- list(sigma_sc = 3)       
    # 	
    if(hasArg(settings)){
    	  #print(suppliedSettings)
          matching 				<- intersect(names(settings), names(suppliedSettings))
          settings[matching] 	<- suppliedSettings[matching]
    }
	scores1 <- rnorm(n, mean = xi_1_mean, sd = sqrt(ev))
	d 		<- scaleDomain2(pnorm(scores1, mean, mean = xi_1_mean, sd = sqrt(ev) ))
	# hist(d)
	#d 		<- scaleDomain( x = 0.75 * scores1 + rnorm(n, sd=sqrt(suppliedSettings$sigma_sc)) ) # 33.5 C: 

	x       <- data.frame(scores1, d)
	colnames(x) <- c("x", "d") # cor(x) cov2cor(sigma_dec)
	return(sim_data 	= x) # plot(d)
}
