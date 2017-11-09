calcFTC <- function(fds, comp_dom, basis_seq , base = "fourier", maxBasisLength = 25, basisChoice = "Med", alpha = 0.05, B = 1000, f_stat = f_stat2, derFds = FALSE){
	### Performs the estimation of the sample
	### Input:
	# - fds: a functional data sample (p x n)
    # - comp_dom: the complete domain (p)
	# - basis_seq: which number of basis length should be checked
    # - base: the kind of basis
	# - maxBasisLength: the maximum allowed basis length for all functions
    # - basisChoice: the statistic applied onto the single basis choices
	# - alpha: significance level for Romano Wolf
    # - B: number of bootstrap replications
	# - f_stat: the regression to obtain statistics applied onto d vs. scores
	# - derFds: derivatives of the functional objects
	### Return:
	# 
	
	## What is the common domain for the fds? fds = simFD$X_ad_m
	observationList <- apply(fds, 2, function(col) which(!is.na(col)))
	small_dom_index	<- min(maxObsIndex <- sapply(observationList, max))
	
	d 				<- comp_dom[maxObsIndex]
	
	small_dom 		<- comp_dom[1:small_dom_index]
	
	krausMean      	<- rowMeans(fds, na.rm = TRUE)
	
	if(base == "fourier"){
		if(!all(basis_seq %% 2 == 1)){
			basis_seq 	<- seq(3, max(basis_seq), 2)
			print("Basis Sequence changed!")
		}
	}
	
	if( !(basisChoice %in% c("Med", "Max")) ){
		stop("Select an overall basis criterion!")
	}
	
	### Choose a number of basis elements b
	bicMat 		<- selBasisLength_sd(X = fds, small_dom, comp_dom, base, basis_seq)
	
	### Select an overall basis length b via Median or Max Basis Length   
	selBasis  	<- get(paste0("calc",basisChoice,"BasisDim"))(bicMat, basis_seq) 		 	
						
	### Restrict maximal basis to 25
	if(selBasis > maxBasisLength){
		selBasis 		<- maxBasisLength
	}
	
	### Take evaluated psi_i and "der" derivative of psi_i of chosen fourier basis of length selBasis base <- "fourier"
	basis_fd_chosen <- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = selBasis)
	# plot(basis_fd_chosen)
	Psi_i			<- eval.basis( evalarg = comp_dom, basisobj = basis_fd_chosen) 
	Psi_i_der		<- eval.basis( evalarg = comp_dom, basisobj = basis_fd_chosen, Lfdobj = 1) 
	
	#matplot(x=comp_dom, y= Psi_i_der[, 75:101], type = "l", col = "darkgrey")
	
	### Estimate coefficients/scores in that chosen basis psi_i, i = 1, ..., selBasis
	coefs <- apply(fds, 2, function(X_i){ # X_i <- fds[, 1]
		#print(which(X_i[1] == fds[1, ]))
		dom_ind 		<- which(!is.na(X_i))
		y               <- as.numeric(na.omit(X_i)) # length(y)
		X_mat           <- Psi_i[dom_ind, ]
		# print(length(dom_ind))
		# In X there is a constant because of Fourier first basis
		return(coefs    <- lm(y ~ -1 + X_mat)$coefficients)
		#plot(y, type = "l")
		#lines(X_mat %*% coefs, col = "red")
	})
	
	
	### Estimate the FTC Estimator like defined in the paper
	ftcMean 		<- calcFTCEstimator(fds, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom, derFds)	
	ftcMean_prime 	<- ftcMean$firstDer
	ftcMean			<- ftcMean$mean
	
    ### Estimate the FTC Estimator with extrapolation which is not in the paper
	#ftcMean2 		<- calcFTCEstimator2(sample = fds, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom)	

	# pdf(paste0(base, "_", dgp, "_MeanAndObsAndMeans.pdf"))
	# matplot(x = comp_dom, y = fds, type = "l", col = addAlpha("darkgrey", alpha = 0.8))
	# lines(x= comp_dom, y = trueMean, type = "l", lwd = 2, col ="black")
	# lines(x= comp_dom, y = kraussMean, type = "l", lwd = 2, col ="blue")
	# lines(x= comp_dom, y = ftcMean, type = "l", lwd = 2, col ="red", lty = 2)
	# dev.off()
	
	### Estimate the FTC Covariance like defined in the paper small_dom = comp_dom
	ftcCov 			<- calcFTCCovariance(fds, ftcMean, ftcMean_prime, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom)	
	#ftcCov2 	<- calcFTCCovariance2(sample = fds, ftcMean = ftcMean, ftcMean_prime = ftcMean_prime, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom)	
	
	krausCov    	<- calcKrausCovariance(fds, ftcMean)
	
	# library(rgl) # dim(ftcCov) str(ftcCov)
	#  persp3d(x = comp_dom, y =comp_dom, z = ftcCov, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = ftcCov, col="skyblue", zlim = c(-5,40))
	# persp3d(x = comp_dom, y =comp_dom, z = ftcCov2, col="skyblue", zlim = c(-5,40))
	# persp3d(x = comp_dom, y =comp_dom, z = trueCov, col="skyblue", zlim = c(-5,40))
	# persp3d(x = comp_dom, y =comp_dom, z = krausCov, col="skyblue", zlim = c(-5,40))
	# persp3d(x = comp_dom, y =comp_dom, z = ftcCov, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = ftcCov2, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = trueCov, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = devCov_est, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = devCov, col="skyblue")
	# persp3d(x = comp_dom, y =comp_dom, z = devCov2, col="skyblue")
	
	### Apply RomanoWolf with F-Statistic
	romWolf <-  RomWolf(X = t(coefs), alpha, B, printout = FALSE, d, f_stat = f_stat2) # plot(X[,1], d)
	
	# matplot(fds, type = "l", col = addAlpha("black", 0.5), ylim = c(-3,3))
	# lines(rowMeans(fds, na.rm = TRUE), col = "red", lwd = 2)
	# lines(mean_sd, col = "blue")
	# abline(h=0, col = "darkgreen", lwd = 2)
	# lines(ftcMean, col = "pink", lwd = 2)
	# mean_sd
	# rowMeans(fds, na.rm = TRUE)
	#}
	retList 	<- list(krausMean = krausMean,
						ftcMean = ftcMean,
						ftcCov = ftcCov,
						krausCov = krausCov,
						bicMat = bicMat,
						coefs = coefs,
						selBasis  = selBasis,
						romWolf = romWolf,
						cor = cor(cbind(d, coefs[1,]))[1,2])
	class(retList) <- append(class(retList),"BID")
	return(retList)

}
