calcKrausCovariance <-
function(fds, fmean){
	########## ########## ########## ########## ########## ##########
	########## Calculation of Covariance matrix according to the definition of Kraus with supplied mean
	# Cov(s,t) = = 1/n sum_i=1^n [(X_i(s) - mu(s))*(X_i(t)-mu(t))]
	
	# sum(sim_data[,"d"] == 1)
	# str(simFD$X_ad_m, max = 1)
    p                   <- length(fmean)

    sample_centered 	<- fds - fmean
    
	# matplot(x = comp_dom, y = sample_centered, type = "l", col = addAlpha("darkgrey", alpha = 0.9))
	
    krausCov			<- matrix(NA, ncol = p, nrow = p)	
   
	# s <- 400;   t <- 500
	# trueCov[s,t]
	# ftcCov[s,t]
    for(s in 1:ncol(krausCov)){ # s = 400 comp_dom[400]
        for(t in 1:ncol(krausCov)){ # t = 400
            X_all_s 			<- sample_centered[ s , ]
            X_all_t 			<- sample_centered[ t , ]
            allNAs 				<- c(which(is.na(X_all_s)),which(is.na(X_all_t)))
            if(length(allNAs) == 0){
                divisor <- length(X_all_s)
                krausCov[s,t] 	<- (X_all_s %*% X_all_t)/divisor
            } else{
                divisor <- length(X_all_s[-allNAs])
                krausCov[s,t] 	<- (X_all_s[-allNAs] %*% X_all_t[-allNAs])/divisor
            }
        }
    }
    return(krausCov = krausCov)
}
