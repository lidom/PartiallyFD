##########################################################################################  
##########################################################################################
## Supplement for   "Partially Observed Functional Data: The Case of Systematically Missing"
## by Dominik Liebl & Stefan Rameseder

rm(list = ls())
library("devtools")
install_github("stefanrameseder/PartiallyFD")
library("PartiallyFD")

##########################################################################################
## Generate Functional Data:
p        		<- 501 # Number of evaluation points per function
comp_dom        	<- seq(from = 0, to = 1,len = p) # Seq. of equidistant evaluation points for complete domain
small_dom		<- comp_dom[ comp_dom <= 0.5 ] # Seq. of equidistant evaluation points for small domain
nbasis_dgp 		<- 5   # Number of basis functions
nbasis_est 		<- 5  # Number of basis functions for estimation process
seq.ev   		<- seq(from = 10, to = 2, len = nbasis_dgp) # Seq. of eigenvalues or variances
sd 				<- sqrt(seq.ev[1])  # standard deviation of normal marginal; in our case the same as the first score
xi_means 		<- c(5,2,0,0,0)
mean_dom		<- 1 # mean of right domain border
prob 			<- 0.9 # binomial probability for normal marginal
basis 			    <- c(BSpline = "bspline", Fourier = "fourier", Monomial = "monomial", Power = "power")
base            <- basis[2] # Just for Fourier
x_var			<- sum(seq.ev)
alpha			<- 0.05 # significance level for bootstraps
B				<- 1000 #number of bootstraps
printout		<- 1
der 			<- 1
base			<- basis[which(basis == "fourier")]
name			<- names(base)

f_stat 			<- PartiallyFD:::f_stat2 # Choose the test statistic function for Romano Wolf

b_choice		<- "Med" # or Max # The choice method for the overall basis selection
maxAllowedBasis	    <- 51 # The number of maximal allowed overall basis length
basis_seq 		<- seq(3,maxAllowedBasis,2) # possible basis lengths
derFds			<- FALSE # Should the estimator be used with the derivatives of the functional data sample or with the basis projection
percentageCompDom   <- 0.02 # how many functions should be forced to be observed over the whole domain
sd_part 		<- 1 # Which part of the d s should be trnaformed to 0.5 or 1 (here 100%)
xi_1_d_cor		<- 0.995 # the correlation in the truncated normal case
d_var 			<- 1 # The variance of d in the truncated normal cases
sigma_sc            <- 3 # the variance of the sigma in the transformation1 case
positions 		<- c(0.625, 0.75, 0.875, 1)
DGP_names		<- c("SD","SC", "ID", "IC") 
ver             <- "FirstSimulation" # Name of dataset to save

# The Supplied Setting list for all DGPs
suppliedSettings <- list(sd_part = sd_part, xi_1_d_cor = xi_1_d_cor, d_var = d_var, sigma_sc = sigma_sc)

replications		<- 1


ns                  <- c(50) # Different n's
#ns                  <- c(50, 150, 250, 500) # Different n's

require(parallel) 		# for parallel spline smoothing

for(n in ns){ # n <- 50
    
    # Version for Saving the .RData
    print(version 			<- paste0(ver,"_p", p, "_n", n, "_R", replications, "_SD", 100*suppliedSettings$sd_part, "_SD", 1000*suppliedSettings$xi_1_d_cor,"_BC_", b_choice,maxAllowedBasis ))
    
    # make the cluster
    cl <- makeCluster(8)
    
    # load required libraries
    clusterEvalQ(cl, {
        library("xtable")
        library("systemfit") 	# install.packages("systemfit") # for SUR estimation
        library("fda")    		# convenient basis functions
        library("scales") 		# colors
        library(copula)			# install.packages("copula")
        library(psych)
        library(vcdExtra)		# install.packages("vcdExtra")
        library(skewt) 			# install.packages("skewt") Skew T Distribution
        library(e1071) 			# for skewness function
        library(MASS) 			# mvrnorm
        library(tmvtnorm)		# install.packages("tmvtnorm")
        library(xtable) 		# latex tables
        library(MBESS)			# install.packages("MBESS")
        library(RColorBrewer) 	# for color definitions
        library(forecast) 		# for arima estimation
        require(MASS) 			# for df
        require(coneproj) 		# cone projections for smoothing parameter lambda
        require(orthogonalsplinebasis) # cone projections for smoothing parameter lambda
        library(PartiallyFD)		# nice corrplots
    })
    
    
    ## export required objects
    clusterExport(cl , varlist = c("p", "n", "comp_dom", "small_dom", "nbasis_dgp", "nbasis_est", "sigma_sc",
                                   "basis_seq", "seq.ev", "sd", "xi_means", "mean_dom", "prob", "b_choice",
                                   "cor", "basis", "x_var", "alpha", "B", "printout", "der", "f_stat", "derFds", "percentageCompDom",
                                   "DGP_names", "replications", "base", "suppliedSettings", "maxAllowedBasis") )
    
    # For each base in Basis
    # for each dgp ind DGPs
    # for all replications
    
    sys.time <- system.time(resDGP <- lapply(DGP_names, function(dgp){ # dgp <- "SC"
        #system.time(resDGP <- lapply(DGP_names, function(dgp){ # dgp <- "ID"
        
        # Define basis once
        basis_fd 		<- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = 5)
        trueMean 		<- PartiallyFD:::calcTrueMean(basis_fd, xi_means, comp_dom) 

        # Calculate true covariance via cov(s,t) = sum_i=1^b lambda_i * psi_i(t)*psi_i(s)
        trueCov 		<- PartiallyFD:::calcTrueCovariance(basis_fd, seq.ev, comp_dom, Lfdobj=0) 
        
        # Calculate true derivative of covariance via cov''(s,t) = sum_i=1^b lambda_i * psi_i'(t)*psi'_i(s)
        trueCovDer 		<- PartiallyFD:::calcTrueCovariance(basis_fd, seq.ev, comp_dom, Lfdobj=1) 
        
        # For all replications
        set.seed(10)
        resRep <- parLapply(cl, 1:replications, function(rep){ # rep = 1
            
            ### Simulate for n cruves xi_i1 and d_i data, where xi_1_mean is the first entry of xi_means 
            if(dgp == "IC"){
                sim_data 					<- PartiallyFD:::simDataIC(n, xi_1_mean = xi_means[1], ev = seq.ev[1], suppliedSettings)
            } else if (dgp == "ID"){
                sim_data 					<- PartiallyFD:::simDataID(n, xi_1_mean = xi_means[1], ev = seq.ev[1], suppliedSettings)
            } else if (dgp == "SC"){
                sim_data 					<- PartiallyFD:::simDataSC(n, xi_1_mean = xi_means[1], ev = seq.ev[1], suppliedSettings)
            } else if (dgp == "SD"){
                sim_data 					<- PartiallyFD:::simDataSD(n, xi_1_mean = xi_means[1], ev = seq.ev[1], suppliedSettings)
            } else {
                print("Could not find simulation function")
            }

            ### If the maximum value of d != 1, define percentageCompDom% of the sample as observed on the full domain
            sim_data[rev(order(sim_data[ ,"d"]))[1:round( percentageCompDom * n ,0)] ,"d"] <- max(comp_dom)
            
            ### Construct corresponding functional data set of n curves observed at p points 
            simFD 						<- PartiallyFD:::simMeanBias(	sim_data = sim_data, nbasis_dgp = nbasis_dgp, p =p, n = n, seq.ev = seq.ev, 
                                        comp_dom = comp_dom, small_dom = small_dom, basis = basis_fd, xi_means = xi_means)
            ### Calculate the Krauss Mean 
            krausMean      <- rowMeans(simFD$X_ad_m, na.rm = TRUE)

            ### Choose a number of basis elements b by Minimizing the BIC over all projections onto basis_seq
            bicMat 		<- PartiallyFD:::selBasisLength(X = simFD$X_ad_m, comp_dom, base = "fourier", basis_seq)
            
            ### Select an overall basis length b via Median or Max Basis Length   
            selBasis  	<- PartiallyFD:::calcMedBasisDim(bicMat, basis_seq) 		 	
            
            ### Take evaluated psi_i and "der" derivative of psi_i of chosen fourier basis of length selBasis
            basis_fd_chosen <- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = selBasis)
            Psi_i			<- eval.basis( evalarg = comp_dom, basisobj = basis_fd_chosen)
            Psi_i_der		<- eval.basis( evalarg = comp_dom, basisobj = basis_fd_chosen, Lfdobj = der)
            
            ### Estimate coefficients/scores in that chosen basis psi_i, i = 1, ..., selBasis
            coefs <- apply(simFD$X_ad_m, 2, function(X_i){ # X_i <- simFD$X_ad_m[, 1]
                dom_ind 		<- which(!is.na(X_i))
                y               <- as.numeric(na.omit(X_i)) # length(y)
                X_mat           <- Psi_i[dom_ind, ]
                return(coefs    <- lm(y ~ -1 + X_mat)$coefficients)
            })
            
            ### Estimate the FTC Estimator like defined in the paper
            ftcMean 		<- PartiallyFD:::calcFTCEstimator(fds = simFD$X_ad_m, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom, derFds = derFds)
            ftcMean_prime 	<- ftcMean$firstDer
            ftcMean			<- ftcMean$mean
            
            ### Estimate the FTC Covariance like defined in the paper small_dom = comp_dom
            ftcCov 		<-  PartiallyFD:::calcFTCCovariance(fds = simFD$X_ad_m, ftcMean = ftcMean, ftcMean_prime = ftcMean_prime, evalBase = Psi_i, evalDer = Psi_i_der, coefs, small_dom, comp_dom)
            krausCov    <-  PartiallyFD:::calcKrausCovariance(fds = simFD$X_ad_m, fmean = krausMean)
            
            if(any(dgp == c("SD", "SC"))){ # dgp = "IC"
                krausCov_ftc    <- PartiallyFD:::calcKrausCovariance(fds = simFD$X_ad_m, fmean = krausMean)
            } else {
                krausCov_ftc    <- krausCov
            }
            
            ### Apply RomanoWolf with F-Statistic
            romWolf <-  PartiallyFD:::RomWolf(X = t(coefs), alpha, B, printout = FALSE, d = sim_data[, "d"], f_stat = f_stat) # plot(X[,1], d)
            
            return(list(krausMean = krausMean,
                        ftcMean = ftcMean,
                        ftcCov = ftcCov,
                        krausCov = krausCov,
                        krausCov_ftc = krausCov_ftc, 
                        selBasis  = selBasis,
                        romWolf = romWolf,
                        cor = cor(sim_data)))
        })	
        return(resRep)
    }))
    
    names(resDGP) <- DGP_names
    stopCluster(cl)	
    print(file <- paste0("PartiallyFD_sim_", version,".RData"))
    print(sys.time)
    save(resDGP, file = paste0("PartiallyFD_sim_", version,".RData"))
}
