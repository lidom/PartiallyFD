
simMeanBias 	<- function(sim_data, nbasis_dgp = 5, p, n, seq.ev,#sim_data= sim_data4
							basis = base, comp_dom, small_dom, xi_means){
	
	# simMeanBias in the continous case
	# sim_data = sim_data_c
	if(any(sim_data[ ,2] < 0)) stop("Negative Domain Value")
	
	## Simulating scores
	scores_m 		<- matrix(NA, nrow = nbasis_dgp, ncol = n) # for each function as column there are nbasis elements for functions

	## First Scores are depending on domain
	scores_m[1, ] 	<- sim_data[ ,1] # in sim_data, the normal marginal have sd = sqrt(seq.ev[b]) 
	
	for(b in 2:nbasis_dgp){					
		scores_m[b, ] 	<- rnorm(n, mean = xi_means[b], sd = sqrt(seq.ev[b]))  
	}
	
	basis		<- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = nbasis_dgp)
	# matplot(eval.basis(evalarg = small_dom, basisobj = basis_est), x = small_dom)
		
	## Save function values on all domains with NAs depending on sim_data[,2]
	X_ad_m 		<- matrix(NA, nrow = p, ncol = n) # for each function as column there are nbasis elements for functions

	for( i in 1:n ){ # i = 101
		# i-th column equals scores times basis as long as domain (comp_c) corresponds
		#print(n)
		domain_i	<- which.max( sim_data[i,"d"] <= comp_dom )  # comp_dom[57] # theoretisch minus 1
		X_ad_m[ 1:(domain_i) ,i] <- eval.basis(evalarg = comp_dom, basisobj = basis)[ 1:domain_i,] %*% scores_m[ ,i]
	}
	
	X_emp_mean		<- rowMeans(X_ad_m, na.rm = TRUE)
	
	dist			<- round(sum(X_emp_mean^2), 2) / n
	
	rownames(scores_m)<- paste0("score",1:nbasis_dgp)
	
	return(list(X_ad_m  = X_ad_m, X_emp_mean = X_emp_mean, basis = basis, 
				scores_m = scores_m, nbasis = nbasis_dgp, p = p, n = n, comp_dom = comp_dom, dist = dist, scores = scores_m))
} 
