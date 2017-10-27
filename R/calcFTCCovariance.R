calcFTCCovariance <-
function(fds, ftcMean, ftcMean_prime, evalBase, evalDer, coefs, small_dom, comp_dom){
	########## ########## ########## ########## ########## ##########
	########## Calculation of Covariance matrix according to the definition in the BID paper
	########## Checks can be done with bias_sim_covarianceEstimation.r
	# Cov(s,t) = = 1/n sum_i=1^n [(X_i(s) - mu(s))*(X_i(t)-mu(t))]
	# estimated via
	# cov_hat(s,t) = 
	# Part I: iff s,t in [0, d_min]
	# 	 	  =	cov_kra(s,t) 
	# Part II: iff s in [0, d_min] and t in (d_min,1]
	# a		  =	int_d_min^t rho(X(s), X'(z))  dz 
	# b       + cov_kra(s,d_min) 
	# Part III:  iff s,t in (d_min,1]
	# a		  =	int_d_min^s int_d_min^t rho(X'(z_2),  X'(z_1)) dz_2 dz_1 
	# b		  + int_d_min^s 			rho(X'(z_1),  X(d_min))     dz_1
	# c		  + int_d_min^t 			rho(X(d_min), X'(z_2))      dz_2
	# d		  + 2*cov_kra(d_min,d_min) 
	# The dimensions of the matrices are with d_1 = length [0, d_min], d_2 = (d_min, 1]
	# Type I:  d_1 x d_1
	# Type II: d_1 x d_2
	# Type III: d_2 x d_2
 	
	
	########## Input:
	# - ftcMean: estimated mean function (px1)
	# - evalBase: evaluated basis (p x b_hat)
	# - evalDer: evaluated derivative of basis (p x b_hat)
	# - coefs: estimated scores per function (n x b_hat)
	# - small_dom: 
	# - comp_dom: 
	########## Output:
	# - cov_fin: the FTC Covariance Matrix of size p x p	
		
	# Number of functions
	n 					<- dim(fds)[2]
	# Number of observation points
	p 					<- length(comp_dom) 
	
	# Specify d_min
	d_min_index			<- which.max(small_dom) # d_1
	
	out_dom 			<- comp_dom[-(1:d_min_index)]
	
	########## ########## ########## ########## ########## ##########
	##########   Part I  ##########
	########## ########## ########## ########## ########## ##########
	########## [0, d_min] 
	# From 0 to d_min where all data is observed; covariance on small domain
	sample_centered 	<- fds[1:d_min_index, ] - ftcMean[1:d_min_index] 
	cov_I				<- (sample_centered %*% t(sample_centered)) / n # dim(sample_centered)
	

	
	# trueMean[d_min_index]; ftcMean[d_min_index]  
	# cov_I[d_min_index,d_min_index]-trueCov[d_min_index,d_min_index]
	# cov_I
	# persp3d(x = small_dom, y = small_dom, z = cov_I-trueCov[1:d_min_index,1:d_min_index], col="skyblue")	
	# library(rgl)
	# persp3d(x = small_dom, y = small_dom, z = cov_I, col="skyblue")	
	# persp3d(x = small_dom, y = small_dom, z = trueCov[1:length(small_dom),1:length(small_dom)], col="skyblue")	
	
	########## ########## ########## ########## ########## ##########
	##########  Part III  ##########
	########## ########## ########## ########## ########## ##########
	########## (d_min,1] ##########
	## From d_min to 1, where some are arbitrary observed; covariance on outer domain

	# Take the first derivatives of the sample
	X_prime 			<- evalDer %*% coefs
	
	# Connect the NAs of the sample with ones in the derivatives
	sample_der 			<- sapply(1:n, function(i){ # i = 10
							# (evalDer %*% coefs)[ , i]
							if(sum(is.na(fds[,i])) != 0){
								v <- c(X_prime[!is.na(fds[,i]), i], rep(NA, times = sum(is.na(fds[,i]))))
								return(v)
							} else {
								return(X_prime[!is.na(fds[,i]), i])
							}
						})
	
	# Center the data, i.e. get  X'_i(t_j) - X'_FTC(t_j) fo for i = 1, ..., n, j = 1, ..., p 
	# dim(sample_der_centered) p x n
	sample_der_centered <- sample_der - matrix(ftcMean_prime, ncol = dim(sample_der)[2], nrow = dim(sample_der)[1]) # dim(sample_der_centered)


	########## 		 Part III a      ##########
	########## (d_min,1] x (d_min,1] ########## 
	# container for derivative of the covariance matrix on the outer domain
	# Cov(s,t) = = 1/n sum_i=1^n [(X_i(s) - mu(s))*(X_i(t)-mu(t))]
	p_outerDom 			<- p - d_min_index
	cov_od_der			<- matrix(NA, ncol = p_outerDom, nrow = p_outerDom)	
	
	for(s in 1:ncol(cov_od_der)){ # s = 1
	 for(t in 1:ncol(cov_od_der)){ # t = 1
		X_all_s 			<- sample_der_centered[ (s+d_min_index) , ]
		X_all_t 			<- sample_der_centered[ (t+d_min_index) , ]
		allNAs 				<- c(which(is.na(X_all_s)),which(is.na(X_all_t)))
		
		if(length(allNAs) == 0){
			divisor <- length(X_all_s)
			cov_od_der[s,t] 	<- (X_all_s %*% X_all_t)/divisor
		} else{
			divisor <- length(X_all_s[-allNAs])
			cov_od_der[s,t] 	<- (X_all_s[-allNAs] %*% X_all_t[-allNAs])/divisor
		}
	 }
	}

	
	# Compare rho'' with rho_hat'' on (d_min,1]^2
	# persp3d(x = out_dom, y =out_dom, z = cov_od_der, col="skyblue")# isSymmetric(cov_hat_III)
	# persp3d(x = out_dom, y =out_dom, z = trueCovDer[(d_min_index+1):p,(d_min_index+1):p], col="skyblue")# isSymmetric(cov_hat_III)
	
	# Integration with respect to the columns
	# cov_od_der_colintegrated =  int rho(X'(t), X'(s)  dt + 
	cov_od_der_colintegrated 	<- apply(cov_od_der, 2, function(col){
										integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = col)
								})
	
	
	# Integration with respect to the rows
	# cov_od = int int rho(X'(t), X'(s) ds dt + 
	cov_hat_IIIa				<- apply(cov_od_der_colintegrated, 1, function(row){
										integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = row)
								})							

	# persp3d(x = 1:250, y = 1:250, z = cov_hat_IIIa, col="skyblue")# isSymmetric(cov_hat_IIIa)
	
	############################################
	########## 	Check with true one:  ##########
	# cov_od_der_colintegrated_true<- apply(trueCovDer[(d_min_index+1):p,(d_min_index+1):p], 2, function(col){
										# integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = col)
								# })
	# cov_hat_IIIa_true			<- apply(cov_od_der_colintegrated_true, 1, function(row){
										# integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = row)
								# })		
	
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = cov_hat_IIIa_true, col="skyblue")# dim(cov_hat_IIIa)
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = cov_hat_IIIa-cov_hat_IIIa_true, col="skyblue")# dim(cov_hat_IIIa)
	# sum((cov_hat_IIIa - cov_hat_IIIa_true))/(250^2)
	############################################
	
	########## 		 Part III b      ##########
	########## int_d_min^s 			rho(X'(z_1),  X(d_min))     dz_1
	# container for mixture covariance matrix rho_hat(X'(z), X(d_min)) with z in in (d_min,1]
	# rho_hat(X'(z), X(d_min)) = sum_i=1^n (X'_i(z) - mu'_i(z))(
	cov_od_mix_der			<-  matrix(NA, ncol = p_outerDom, nrow = 1)	
	
	# X(d_min_index) - mu_ftc(d_min_index)
	X_all_t 			<- sample_centered[ d_min_index , ] 
	# simFD$X_ad_m[d_min_index ,] - ftcMean[d_min_index]
	# dim(sample_der_centered)
	for(s in 1:ncol(cov_od_mix_der)){ # s = 1
		
		# X'(z)-mu'(z)
		X_all_s				<- sample_der_centered[ (s+d_min_index) , ] # length(X_all_s) 
		allNAs 				<- c(which(is.na(X_all_s)),which(is.na(X_all_t)))
		# length(X_all_s[-allNAs])
		if(length(allNAs) == 0){
			divisor <- length(X_all_s)
			cov_od_mix_der[1,s] 	<- (X_all_s %*% X_all_t)/divisor
		} else{
			divisor <- length(X_all_s[-allNAs])
			cov_od_mix_der[1,s] 	<- (X_all_s[-allNAs] %*% X_all_t[-allNAs])/divisor
		}
	} 
	
	############################################
	########## 	Check with true one:  ##########
	# trueMixCov <- calcTrueMixedCovariance(basis_fd, seq.ev, comp_dom); trueMixCov[d_min_index, ] == trueMixCov[ ,d_min_index]
	# plot( x= out_dom, y= c(cov_od_mix_der), type = "l")
	# lines(x = comp_dom, y = trueMixCov, type = "l", col = "red")
	
	cov_hat_IIIb_border 			<- integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = cov_od_mix_der)
	cov_hat_IIIb 					<- matrix(cov_hat_IIIb_border, ncol = p_outerDom, nrow = p_outerDom)
	# plot( x= out_dom, y= c(cov_od_mix_der), type = "l")
	# lines(x = comp_dom, y = trueMixCov, type = "l", col = "red")
	############################################
	cov_hat_IIIc 				<- t(cov_hat_IIIb)
		
	
	########## 		 Part III d      ##########
	########## cov(d_min, d_min)
	cov_hat_IIId			<- cov_I[d_min_index,d_min_index]
	
	# Check with true:
	# trueCov[d_min_index,d_min_index]
	
	########## 	 Part III total     ##########
	cov_III					<- cov_hat_IIIa + cov_hat_IIIb + cov_hat_IIIc + cov_hat_IIId
	# sum((cov_III - trueCov[(d_min_index+1):p,(d_min_index+1):p])^2)/(250^2)
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = round(abs(cov_III - trueCov[(d_min_index+1):p,(d_min_index+1):p]),1), col="skyblue", zlim = c(0, 8))
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = cov_III, col="skyblue", zlim = c(-5, 35))# isSymmetric(cov_hat_III)
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = trueCov[(d_min_index+1):p,(d_min_index+1):p], zlim = c(-5, 35), col="green")# isSymmetric(cov_hat_III)
	# cov_hat_IIIa[1,250]
	# cov_hat_IIIb[1,250]
	# cov_hat_IIIc[1,250]
	# cov_hat_IIId
	# library(rgl)
	# persp3d(x = (d_min_index+1):p, y =(d_min_index+1):p, z = cov_III, col="skyblue")# isSymmetric(cov_hat_III)
	# persp3d(x = (d_min_index+1):p, y = (d_min_index+1):p, z = trueCov[(d_min_index+1):p,(d_min_index+1):p] , col="skyblue")# isSymmetric(cov_hat_III)
	
	########## ########## ########## ########## ########## ##########
	##########   Part II  ##########
	########## ########## ########## ########## ########## ##########
	########## 		 Part II a      ##########	
	########## [0,d_min] x (d_min,1] ########## 
	# container for mixture covariance matrix rho(X(s), X'(z)) with s in [0, d_min] and t in (d_min,1]
	# in total a d_1 x d_2 matrix
	# covariance for Mixture Domain where the one on the outer is in the derivative
	cov_md_der			<- matrix(NA, ncol = p_outerDom   , nrow = d_min_index)	
	
	for(s in 1:nrow(cov_md_der)){ # s = 16
	 for(t in 1:ncol(cov_md_der)){ # t = 10
		X_all_s 			<- sample_centered[ s , ]
		X_all_t 			<- sample_der_centered[ (t+d_min_index) , ]
		allNAs 				<- c(which(is.na(X_all_s)),which(is.na(X_all_t)))
		# length(X_all_s[-allNAs])
		if(length(allNAs) == 0){
			divisor <- length(X_all_s)
			cov_md_der[s,t] 	<- (X_all_s %*% X_all_t)/divisor
		} else{
			divisor <- length(X_all_s[-allNAs])
			cov_md_der[s,t] 	<- (X_all_s[-allNAs] %*% X_all_t[-allNAs])/divisor
		}
	 }
	}
	# dim(cov_md_der) # rho(X'(t),X(s)) on (d_min,1] x [0, d_min]

	
	# Integration with respect to t in (d_min,1], i.e. wrt. to rows
	# i.e. cov_hat_IIa	 = int_d_min^t rho(X(s), X'(z))  dz 
	
	cov_hat_IIa				<- t(apply(cov_md_der, 1, function(col){ # row = 
										integrate_xy_left(x = comp_dom[ (d_min_index+1):p ], y = col)
								})) # dim(cov_hat_IIa); length((d_min_index+1):p)
	
	########## 		 Part II b      ##########
	########## rho(s, d_min) for s in [0,d_min] // cov_kra(s,d_min)
	# in total a d_1 x d_2 matrix
	cov_hat_IIb			<- matrix(cov_I[1:d_min_index , d_min_index], ncol = p_outerDom, nrow = d_min_index )
	
	
	########## 	 Part II total      ##########
	cov_II 				<- cov_hat_IIa + cov_hat_IIb # dim(cov_hat_IIa); dim(cov_hat_IIb)
	
	
	########## ########## ########## ########## ########## ##########
	########## 	 Connect Cov Matrix to get a p x p Matrix  ##########
	cov_fin 			<- matrix(NA, ncol = p, nrow = p)
	
	cov_fin[1:d_min_index, 1:d_min_index] 			<- cov_I
	cov_fin[(d_min_index+1):p, (d_min_index+1):p] 	<- cov_III
	cov_fin[1:d_min_index, (d_min_index+1):p] 		<- cov_II
	cov_fin[(d_min_index+1):p,1:d_min_index] 		<- t(cov_II)
	
	# cov_fin[(d_min_index-2):(d_min_index+2),(d_min_index-2):(d_min_index+2)] 
	# cov_fin[,251:252] 
	# cov_fin[501,251:252] # difference of 16
	
	# isSymmetric(cov_fin)
	# library(rgl)
	# persp3d(x = comp_dom, y = comp_dom, z = cov_fin, col="skyblue", zlim=c(-15,45))# isSymmetric(cov_hat_III)
	# persp3d(x = comp_dom, y = comp_dom, z = krausCov, col="green", zlim=c(-15,45))# isSymmetric(cov_hat_III)
	# rgl.snapshot( "estimatedCov.png", fmt = "png", top = TRUE )
	# persp3d(x = comp_dom, y = comp_dom, z = trueCov, col="red", zlim=c(-15,45))# isSymmetric(cov_hat_III)
	# rgl.snapshot( "trueCov.png", fmt = "png", top = TRUE )
	# persp3d(x = comp_dom, y = comp_dom, z = cov_fin- trueCov, col="skyblue")# isSymmetric(cov_hat_III)
	# rgl.snapshot( "estMinusTrueCov.png", fmt = "png", top = TRUE )

	return(cov_fin = cov_fin)
	
}
