RomWolf <-
function(X, alpha, B, printout, d, f_stat){
    k       	<- dim(X)[2] # number of hypotheses
    n       	<- dim(X)[1] # number of observations
    
    # compute t statistics for each hypothesis d ~ beta score j and order them 
    #t_statv 	<- t(apply(X, 2 , f_stat, d = d )) 
    t_statv 	<- f_stat(X, d, k) # plot(X[,1], d) f_stat = f_stat1
	
    # sometimes there are NAs because of linear dependencies, then one have to adjust
    # the number of k
    # k           <- length(t_statv)
    t_stat_rk 	<- rank(-t_statv)  # ? rank(c(1,2,5,-10))
	
	# if (printout) {
	# 	cat("t statistics ", t_statv, "\n")
	# 	cat("Ranks of t statistics ", t_stat_rk, "\n")  
	# }
    
	# compute Bootstrap p values
    t_bootm 	<- matrix(NA, ncol=k, nrow=B)
    X_mean  	<- colMeans(X) 
	
	d_star 		<- matrix(summary(lm(d ~ 1 + ., data = data.frame(d, X)))$resid, ncol = 1)
    
	Km      	<- matrix(1:k,ncol=k) %x% matrix(1,nrow=k,ncol=1) # matrix to store remaining hypotheses
	
	# f_stat <- f_stat2
	#set.seed(1)
	#system.time(
	for (b in 1:B){ # b = 1
	    #print(b)
        bootcol     <- sample.int(n, replace = TRUE)
		d_boot		<- d_star[bootcol,]
        X_boot      <- matrix(X[bootcol,], nrow = n)  - matrix(X_mean, ncol = k) %x% matrix(1,ncol = 1,nrow=n)		# dim(X_boot)
		t_bootm[b,] <- f_stat( X = X_boot, d = d_boot, k)
    }   
	
	
	
    for (j in (1:k)){ # j <- 1 # for all hypotheses
        K       	<- na.omit(Km[j,]) # K th row
		
        # compute Bootstrap p value
        B_here 		<- B-1
        
		if (length(K)>1){
            max_t_stat_sample <- max(t_statv[K]) %x% matrix(1,nrow=B_here,ncol=1)
            max_t_stat_boot   <- apply(t_bootm[(1:B_here),K],1,max) 
        } else {
            max_t_stat_sample <- max(t_statv[K]) %x% matrix(1,nrow=B_here,ncol=1)
            max_t_stat_boot   <- t_bootm[(1:B_here)]
        }
            
        pval 		<- ( 1 + sum(max_t_stat_boot >= max_t_stat_sample) ) / B
        #pval <- mean( (max_t_stat_boot >= max_t_stat_sample) )                
        
        #if (printout) cat("p value is ", pval, "\n")
        if (pval < alpha){
            H_reject    		<- which(t_stat_rk==j)
            if (printout) cat("H_reject ", H_reject, "\n")
            Km[(j:k),H_reject]  <- NA
            if (printout) cat("Hypothesis", H_reject, " is rejected.\n")
        } else {
            if (printout) cat("Intersection hypothesis not rejected.\n")
            break
        }
    }
    ent <- (is.na(Km[k,])) # puts 1 if hypothesis j is rejected
    #return(list(ent=ent,Km=Km,pval=pval, resid = d_star))
    return(list(ent=ent,pval=pval))
}
