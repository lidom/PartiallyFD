##########################################################################################  
##########################################################################################
## Supplement for "Partially Observed Functional Data: The Case of Systematically Missing"
## by Dominik Liebl & Stefan Rameseder						            
## Analyse the results from the predecessing simulation


# Need parameters for version since all files saved before will be loaded here
version     <- paste0(ver,"_p", p, "_n", n, "_R", replications, "_SD", 100*suppliedSettings$sd_part, "_SD", 1000*suppliedSettings$xi_1_d_cor,"_BC_", b_choice,maxAllowedBasis ))

maxAllowedBasis	    <- 51
replications	    <- 500
ns                  <- c(50,150,250,500)


require(xtable)
library(xtable)
##########################################################################################
## Get it for all n's together into one list
allResults <- lapply(ns, function(x){ #
    
    print(version 			<- "FirstSimulation_p501_n50_R1_SD100_SD995_BC_Med31")
    load(paste0("PartiallyFD_sim_", version,".RData"))
    DGP_names   <- names(resDGP)
    cor         <- lapply(resDGP, function(dgp) sapply(dgp, function(rep) rep$cor[1,2]))

    (corSum <- lapply(cor, summary))
    
    ##########################################################################################
    ### a) Mean Estimation
    ## True Mean
    # Define basis once
    basis_fd 		<- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = nbasis_dgp)
    
    # Calculate true mean
    trueMean 		<- calcTrueMean(basis_fd, xi_means, comp_dom) 
    
    
    ## FTC Estimator
    ftcMeans    	<- lapply(resDGP, function(dgp) lapply(dgp, function(rep) rep$ftcMean))
    
    ## Bias: Calculate the Mean of all Means
    # Bias[\hat{\mu}_{n}(t)] = 
    #                        = E_{MC}[\hat{\mu}_{n}(t)  ] - \mu(t) 
    #                        = (1/500)\sum_{r=1}^{500}[\hat{\mu}_{n,r}(t)] - \mu(t)
    ftcMeansMeans   <- lapply(ftcMeans, function(dgp){
                           return(rowMeans(matrix(unlist(dgp), nrow = 501) ))
    } )
    ftcMeansBias    <- lapply(ftcMeansMeans, function(rep) rep - trueMean)
    SqIntMeanBias   <- lapply(ftcMeansBias, function(rep) sum((rep)^2)/p)
    
    ## Calculate the Variance of the Means
    
    #Var[\hat{\mu}_{n}(t)] =
    #                      = E_{MC}[ (\hat{\mu}_{n}(t) - E_{MC}[\hat{\mu}_{n}(t)  ] )^2 ]
    #                      = (1/500)\sum_{r=1}^{500}[ (\hat{\mu}_{n,r}(t) - (1/500)\sum_{l=1}^{500}[\hat{\mu}_{n,l}(t)] )^2 ] 
    VarOfMeans   <- vector("list", length = length(DGP_names))
    names(VarOfMeans) <- DGP_names
    for(dgp in DGP_names){ # dgp <- DGP_names[[1]]; rep <- ftcMeans[[1]][[1]]
        VarOfMeans[[dgp]] <- lapply(ftcMeans[[dgp]], function(rep) ( rep - ftcMeansMeans[[dgp]] )^2 )
    }
    ftcVarMeans   <- lapply(VarOfMeans, function(dgp){
        return(rowMeans(matrix(unlist(dgp), nrow = p) ))
    } )
    IntftcVarMeans   <- lapply(ftcVarMeans, function(rep) sum(rep)/p)
    
    ## MSE Median and Mean
    sqDevs 			<- lapply(ftcMeans, function(dgp) sapply(dgp, function(rep) sum((rep-trueMean)^2)/length(comp_dom)))
    ftcSqDevsSum 	<- lapply(sqDevs, function(dgp) summary(dgp))
    ftcMeansAndMed  <- t(sapply(ftcSqDevsSum, function(x) x[c("Mean", "Median")]))
    
    ## Krauss Estimator
    krausMeans 		<- lapply(resDGP, function(dgp) lapply(dgp, function(rep) rep$krausMean))
    
    ## Bias: Calculate the Mean of all Means
    krausMeansMeans   <- lapply(krausMeans, function(dgp){
        return(rowMeans(matrix(unlist(dgp), nrow = p) ))
    } )
    krausMeansBias    <- lapply(krausMeansMeans, function(rep) rep - trueMean)
    krausSqIntMeanBias   <- lapply(krausMeansBias, function(rep) sum((rep)^2)/p)
    
    ## Calculate the Variance of the Means
    VarOfMeans   <- vector("list", length = 4)
    names(VarOfMeans) <- DGP_names
    for(dgp in DGP_names){ # dgp <- DGP_names[[4]]; rep <- ftcMeans[[1]][[1]]
        VarOfMeans[[dgp]] <- lapply(krausMeans[[dgp]], function(rep) ( rep - krausMeansMeans[[dgp]] )^2 )
    }
    krausVarMeans   <- lapply(VarOfMeans, function(dgp){
        return(rowMeans(matrix(unlist(dgp), nrow = p)) )
    } )
    IntkrausVarMeans   <- lapply(krausVarMeans, function(rep) sum(rep)/p)

    krausptwDevs    <- lapply(ftcMeans, function(dgp) sapply(dgp, function(rep) rep-trueMean))
    krausSqDevs 	<- lapply(krausMeans, function(dgp) sapply(dgp, function(rep) sum((rep-trueMean)^2)/length(comp_dom)))
    krausSqDevsSum 	<- lapply(krausSqDevs, function(dgp) summary(dgp))
    krausMeansAndMed<- t(sapply(krausSqDevsSum, function(x) x[c("Mean", "Median")]))
    
    ## Bring both together
    EST_names   	<- c("ftc", "kraus")
    distMatrix 		<- matrix(NA, ncol = 2*length(DGP_names), nrow = 2)
    colnames(distMatrix) <- paste(rep(c("Mean", "Median"), length(DGP_names)), rep(DGP_names, each =2))
    rownames(distMatrix) <- EST_names
    
    for(est in EST_names){ # est <- ftc
        for(stat in c("Mean", "Median")){ # stat <- "Mean"
            for(dgp in DGP_names){ # dgp <- "SD"
                #print(paste0(est, stat, dgp))
                distMatrix[est,paste(stat,dgp)] <- get(paste0(est,"MeansAndMed"))[dgp,stat]
            }
        }
    }
    
    (distMatrixLat <- round(distMatrix,2))
    
    ## FTC Estimator
    ftcMeans    	<- lapply(resDGP, function(dgp) lapply(dgp, function(rep) rep$ftcMean))
    sqDevs 			<- lapply(ftcMeans, function(dgp) sapply(dgp, function(rep) sum((rep[-(1:250)]-trueMean[-(1:250)])^2)/length(comp_dom[-(1:250)])))
    ftcSqDevsSum 	<- lapply(sqDevs, function(dgp) summary(dgp))
    ftcMeansAndMed  <- t(sapply(ftcSqDevsSum, function(x) x[c("Mean", "Median")]))
    
    ## Krauss Estimator length(-(251:501))
    krausMeans 		<- lapply(resDGP, function(dgp) lapply(dgp, function(rep) rep$krausMean))
    krausSqDevs 	<- lapply(krausMeans, function(dgp) sapply(dgp, function(rep) sum((rep[-(1:250)]-trueMean[-(1:250)])^2)/length(comp_dom[-(1:250)])))
    krausSqDevsSum 	<- lapply(krausSqDevs, function(dgp) summary(dgp))
    krausMeansAndMed<- t(sapply(krausSqDevsSum, function(x) x[c("Mean", "Median")]))
    
    
    distMatrixsd 		<- matrix(NA, ncol = 2*length(DGP_names), nrow = 2)
    colnames(distMatrixsd) <- paste(rep(c("Mean", "Median"), length(DGP_names)), rep(DGP_names, each =2))
    rownames(distMatrixsd) <- EST_names
    
    for(est in EST_names){ # est <- ftc
        for(stat in c("Mean", "Median")){ # stat <- "Mean"
            for(dgp in DGP_names){ # dgp <- "SD"
                #print(paste0(est, stat, dgp))
                distMatrixsd[est,paste(stat,dgp)] <- get(paste0(est,"MeansAndMed"))[dgp,stat]
            }
        }
    }
    
    (distMatrixsdLat <- round(distMatrixsd,2))

    ##########################################################################################
    ### b) Chosen Basis Length
    selBasis    <- lapply(resDGP, function(dgp) sapply(dgp, function(rep) rep$selBasis))
    lapply(selBasis, table)
    
    ##########################################################################################
    ### c) Identification via Romano Wolf
    romWolfent    <- lapply(resDGP, function(dgp) lapply(dgp, function(rep) rep$romWolf$ent))
    
    # Define scenario F1, F2, and the Correct One
    # F1 = All are not rejected
    # F2 = There are more than the first rejected
    # C = Only H1 is rejected, the others are not
    
    romWolfDec <- lapply(romWolfent, function(dgp) sapply(dgp, checkFtcHypothesis))
    lapply(romWolfDec, table)
    ##########################################################################################
    ### d)  Analysis of Correlations
    cor    <- lapply(resDGP, function(dgp) sapply(dgp, function(rep) rep$cor[1,2]))
    
    # lapply(resDGP[["SC"]], function(x) x$cor[1,2])
    (corSum <- lapply(cor, summary))
    
    ##########################################################################################
    ### e)  Analysis of Covariances
    resCov 	<- lapply(DGP_names, function(dgp){ # dgp <- "SC"
        # Define basis once
        basis_fd 		<- get(paste0("create.",base,".basis"))(rangeval = range(comp_dom), nbasis = nbasis_dgp)
        
        # Calculate true mean
        trueMean 		<- calcTrueMean(basis_fd, xi_means, comp_dom) 
        
        # Calculate true covariance via cov(s,t) = sum_i=1^b lambda_i * psi_i(t)*psi_i(s)
        trueCov 		<- calcTrueCovariance(basis_fd, seq.ev, comp_dom) 
        
        return(list(trueMean = trueMean, trueCov = trueCov))
    })
    names(resCov) <- DGP_names
    
    covDistMatrix <- matrix(NA, ncol = dim(distMatrixLat)[2], nrow = dim(distMatrixLat)[1])
    dimnames(covDistMatrix) <- dimnames(distMatrix)
 
    for(dgp in DGP_names){ # dgp <- "SD"
        print(paste0(dgp))
        ftcCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$ftcCov) # rep <- ftcCovs[[1]]
        ftcCovDist 		<- sapply(ftcCovs, function(rep) sum((rep - resCov[[dgp]]$trueCov)^2)/p^2)
        covDistMatrix["ftc", paste(c("Mean", "Median"),dgp) ] <- c(mean(ftcCovDist), median(ftcCovDist))
        
        krausCovs    	<- lapply(resDGP[[dgp]], function(rep) rep$krausCov)
        ftcCovDist 		<- sapply(krausCovs, function(rep) sum((rep - resCov[[dgp]]$trueCov)^2)/p^2)
        covDistMatrix["kraus",paste(c("Mean", "Median"),dgp)] <- c(mean(ftcCovDist), median(ftcCovDist))
    }
    
    (covDistMatrixLat <- round(covDistMatrix,2))

    covDistMatrixsd <- matrix(NA, ncol = dim(distMatrixLat)[2], nrow = dim(distMatrixLat)[1])
    dimnames(covDistMatrixsd) <- dimnames(distMatrix)
    
    for(dgp in DGP_names){ # dgp <- "SD"
        print(paste0(dgp))
        ftcCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$ftcCov) # rep <- ftcCovs[[1]]
        ftcCovDist 		<- sapply(ftcCovs, function(rep) sum((rep[-(1:250), -(1:250)] - resCov[[dgp]]$trueCov[-(1:250), -(1:250)])^2)/(p-250)^2)
        covDistMatrixsd["ftc", paste(c("Mean", "Median"),dgp) ] <- c(mean(ftcCovDist), median(ftcCovDist))
        
        krausCovs    	<- lapply(resDGP[[dgp]], function(rep) rep$krausCov)
        ftcCovDist 		<- sapply(krausCovs, function(rep) sum((rep[-(1:250), -(1:250)] - resCov[[dgp]]$trueCov[-(1:250), -(1:250)])^2)/(p-250)^2)
        covDistMatrixsd["kraus",paste(c("Mean", "Median"),dgp)] <- c(mean(ftcCovDist), median(ftcCovDist))
    }
    
    (covDistMatrixsdLat <- round(covDistMatrixsd,2))
    
	covDistMatrix_ftc <- matrix(NA, ncol = dim(distMatrixLat)[2], nrow = dim(distMatrixLat)[1])
    dimnames(covDistMatrix_ftc) <- dimnames(distMatrix)
 
    for(dgp in DGP_names){ # dgp <- "SD"
        print(paste0(dgp))
        ftcCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$ftcCov) # rep <- ftcCovs[[1]]
        ftcCovDist 		<- sapply(ftcCovs, function(rep) sum((rep - resCov[[dgp]]$trueCov)^2)/p^2)
        covDistMatrix_ftc["ftc", paste(c("Mean", "Median"),dgp) ] <- c(mean(ftcCovDist), median(ftcCovDist))
        
        krausCovs    	<- lapply(resDGP[[dgp]], function(rep) rep$krausCov_ftc)
        ftcCovDist 		<- sapply(krausCovs, function(rep) sum((rep - resCov[[dgp]]$trueCov)^2)/p^2)
        covDistMatrix_ftc["kraus",paste(c("Mean", "Median"),dgp)] <- c(mean(ftcCovDist), median(ftcCovDist))
    }
    (covDistMatrixLat_ftc <- round(covDistMatrix_ftc,2))
    
	covDistMatrixsd_ftc <- matrix(NA, ncol = dim(distMatrixLat)[2], nrow = dim(distMatrixLat)[1])
    dimnames(covDistMatrixsd_ftc) <- dimnames(distMatrix)
    
    
    for(dgp in DGP_names){ # dgp <- "SD"
        print(paste0(dgp))
        ftcCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$ftcCov) # rep <- ftcCovs[[1]]
        ftcCovDist 		<- sapply(ftcCovs, function(rep) sum((rep[-(1:250), -(1:250)] - resCov[[dgp]]$trueCov[-(1:250), -(1:250)])^2)/(p-250)^2)
        covDistMatrixsd_ftc["ftc", paste(c("Mean", "Median"),dgp) ] <- c(mean(ftcCovDist), median(ftcCovDist))
        
        krausCovs    	<- lapply(resDGP[[dgp]], function(rep) rep$krausCov_ftc)
        ftcCovDist 		<- sapply(krausCovs, function(rep) sum((rep[-(1:250), -(1:250)] - resCov[[dgp]]$trueCov[-(1:250), -(1:250)])^2)/(p-250)^2)
        covDistMatrixsd_ftc["kraus",paste(c("Mean", "Median"),dgp)] <- c(mean(ftcCovDist), median(ftcCovDist))
    }
    
    (covDistMatrixsdLat_ftc <- round(covDistMatrixsd_ftc,2))

    covDistBias <- matrix(NA, ncol = 2*length(DGP_names), nrow = 2)
    dimnames(covDistBias) <- dimnames(distMatrix)
    for(dgp in DGP_names){ # dgp <- DGP_names[1]
        # FTC
        # a) Bias
        ## Bias: Calculate the Mean of all Means
        # Bias[\hat{\sigma}_{n}(s,t)] = 
        #                        = E_{MC}[\hat{\sigma}_{n}(s,t)  ] - \hat{\sigma}_{n}(s,t) 
        #                        = (1/500)\sum_{r=1}^{500}[\hat{\sigma}_{n,r}(s,t)] - \sigma(s,t)
        
        print(dgp)
        ftcCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$ftcCov) # rep <- ftcCovs[[1]]
        #str(ftcCovs, max = 1); length(ftcCovs) # rep <- ftcCovs[[1]]
        #str(ftcCovs1, max = 1); length(ftcCovs1) # rep <- ftcCovs[[1]]
        
        # Stack matrix to vector
        ftcCovs1        <- lapply(ftcCovs, function(rep) c(rep))
        # Calculate rowMeans over all stacked matrices
        ftcCovs2        <- rowMeans(matrix(unlist(ftcCovs1), nrow = p^2))
        # Back transform to matrix
        ftcCovs3         <- matrix(ftcCovs2, byrow = TRUE, nrow = p)
        
        ftcCovBias <-  sum( (ftcCovs3 - resCov[[dgp]]$trueCov)^2 ) /p^2
        
        # b) Var
        #Var[\hat{\sigma}_{n}(t)] =
        #                      = E_{MC}[ (\hat{\sigma}_{n}(t) - E_{MC}[\hat{\sigma}_{n}(t)  ] )^2 ]
        #                      = (1/500)\sum_{r=1}^{500}[ (\hat{sigma}_{n,r}(t) - (1/500)\sum_{l=1}^{500}[\hat{sigma}_{n,l}(t)] )^2 ] 
        #                      = (1/500)\sum_{r=1}^{500}[ (\hat{\sigma}_{n,r}(t) - ftcCovs3 )^2 ] 
        ftcCovs_Var <- lapply(ftcCovs, function(rep) (rep - ftcCovs3)^2)
        
        ftcCovs_Var1        <- lapply(ftcCovs_Var, function(rep) c(rep))
        # Calculate rowMeans over all stacked matrices
        ftcCovs_Var2        <- rowMeans(matrix(unlist(ftcCovs_Var1), nrow = p^2))
        # Back transform to matrix
        ftcCovs_Var3         <- matrix(ftcCovs_Var2, byrow = TRUE, nrow = p)
        ftcCovVar <- sum(ftcCovs_Var3)/p^2
        #c) Insert into matrix
        covDistBias["ftc",paste(c("Mean", "Median"),dgp)] <- c(ftcCovBias,ftcCovVar)
        
        # Kraus
        # a) Bias
        krausCovs    		<- lapply(resDGP[[dgp]], function(rep) rep$krausCov) # rep <- ftcCovs[[1]]
        
        # Stack matrix to vector
        ftcCovs1        <- lapply(krausCovs, function(rep) c(rep))
        # Calculate rowMeans over all stacked matrices
        ftcCovs2        <- rowMeans(matrix(unlist(ftcCovs1), nrow = p^2))
        # Back transform to matrix
        ftcCovs3         <- matrix(ftcCovs2, byrow = TRUE, nrow = p)
        
        krausCovBias <-  sum( (ftcCovs3 - resCov[[dgp]]$trueCov)^2 ) /p^2
        
        # b) Var
        #Var[\hat{\sigma}_{n}(t)] =
        #                      = E_{MC}[ (\hat{\sigma}_{n}(t) - E_{MC}[\hat{\sigma}_{n}(t)  ] )^2 ]
        #                      = (1/500)\sum_{r=1}^{500}[ (\hat{sigma}_{n,r}(t) - (1/500)\sum_{l=1}^{500}[\hat{sigma}_{n,l}(t)] )^2 ] 
        #                      = (1/500)\sum_{r=1}^{500}[ (\hat{\sigma}_{n,r}(t) - ftcCovs3 )^2 ] 
        krausCovs_Var <- lapply(krausCovs, function(rep) (rep - ftcCovs3)^2)
        
        ftcCovs_Var1        <- lapply(krausCovs_Var, function(rep) c(rep))
        # Calculate rowMeans over all stacked matrices
        ftcCovs_Var2        <- rowMeans(matrix(unlist(ftcCovs_Var1), nrow = p^2))
        # Back transform to matrix
        ftcCovs_Var3         <- matrix(ftcCovs_Var2, byrow = TRUE, nrow = p)
        
        krausCovVar <- sum(ftcCovs_Var3)/p^2
        
        # c) Insert into matrix
        covDistBias["kraus",paste(c("Mean", "Median"),dgp)] <- c(krausCovBias,krausCovVar)
    }
    colnames(covDistBias) <- paste(rep(c("Bias", "Variance"), length(DGP_names)), rep(DGP_names, each =2))
    ##########################################################################################
    ##########################################################################################
    ### Summary
    (distMatrixLat <- round(distMatrix,2))
    corSum
    lapply(selBasis, table)
    lapply(lapply(romWolfDec, table), function(x) round(x/replications*100,1))
    (covDistMatrixLat <- round(covDistMatrix,2))
    
    return(list(distMatrixLat = distMatrixLat, distMatrixsdLat = distMatrixsdLat,
                corSum = corSum, selBasis = selBasis, romWolfDec = romWolfDec, covDistMatrixsdLat = covDistMatrixsdLat,
                covDistMatrixLat = covDistMatrixLat, covDistMatrixsdLat_ftc = covDistMatrixsdLat_ftc,
                ftcBias = SqIntMeanBias,
                krausBias = krausSqIntMeanBias,
                ftcVar = IntftcVarMeans,
                krausVar = IntkrausVarMeans,
                covDistBias = covDistBias,
                covDistMatrixLat_ftc = covDistMatrixLat_ftc))
})

names(allResults ) <-  ns
save(allResults, file = paste0("R_Data/bid_sim_complete.RData"))



##########################################################################################
#### In all allResults now a list of
# n
## dgp
load(file = paste0("R_Data/bid_sim_complete.RData"))

path <- "/Users/stefanrameseder/Dropbox/BID DP/Round_2_CSDA_SI_Due_Apr_14_2018/Manuscript/bid_sim_means.tex"
path <- "C:\\Users\\LocalAdmin\\Dropbox\\BID DP\\Manuscript\\latexChapters\\bid_sim_means.tex"

## Latex Options
table_placement    	<- "!htpb"
hlineAfter			<- seq(2,8, 2)
# Means
distMat 			<- lapply(allResults, function(n) n$distMatrixLat)
(tdistMat			<- transfBidSimObj(distMat))

print(xtable(tdistMat) , type="latex", file=paste0(path), include.rownames = TRUE, 
	  hline.after = hlineAfter,size = "\\setlength{\\tabcolsep}{0.1cm}", 
      table.placement = getOption("xtable.table.placement", paste(table_placement)))

	  
distMatSD 			<- lapply(allResults, function(n) n$distMatrixsdLat)
tdistMatSD			<- transfBidSimObj(distMatSD)

print(xtable(tdistMatSD) , type="latex", file=path, include.rownames = TRUE, 
	  hline.after = hlineAfter,size = "\\setlength{\\tabcolsep}{0.1cm}", 
      table.placement = getOption("xtable.table.placement", paste(table_placement)))

# Biases and Variances

ftcbias <- lapply(allResults, function(n) n$ftcBias)
krausbias <- lapply(allResults, function(n) n$krausBias)
ftcvar <- lapply(allResults, function(n) n$ftcVar)
krausvar <- lapply(allResults, function(n) n$krausVar)
str(ftcBias, max = 2)
biasAndVars		<- matrix(NA, ncol = length(DGP_names)*2, nrow = length(ns)*2)
colnames(biasAndVars) <- paste0(c("bias", "var"),rep(DGP_names, each = 2))
rownames(biasAndVars) <-  paste0(c("ftc", "kraus"),rep(ns, each = 2))
for(n in ns){ # n <- ns[1]
    for(est in c("ftc", "kraus")){ # est <- "ftc"
        for(stat in c("bias", "var")){ # stat <- "bias"
            for(dgp in DGP_names){ # dgp <- DGP_names[1]
                biasAndVars[ paste0(est,as.character(n)) , paste0(stat, dgp)] <- get(paste0(est, stat))[[as.character(n)]][[dgp]]
            }
        }
    }
}
round(biasAndVars, 3)
unlist(ftcbias)
unlist(ftcvar)
# Bring all together



# Covariances
covDistMatrix 		<- lapply(allResults, function(n) n$covDistBias)
(tcovDistMatrix		<- round(transfBidSimObj(covDistMatrix),1))

print(xtable(tcovDistMatrix) , type="latex", file=path, include.rownames = TRUE, 
	  hline.after = hlineAfter,size = "\\setlength{\\tabcolsep}{0.1cm}", 
      table.placement = getOption("xtable.table.placement", paste(table_placement)))

covDistMatrixSD 	<- lapply(allResults, function(n) n$covDistMatrixsdLat)
tcovDistMatrixSD	<- transfBidSimObj(covDistMatrixSD)

print(xtable(tcovDistMatrixSD) , type="latex", file=paste0("C:\\Users\\LocalAdmin\\Dropbox\\BID DP\\Manuscript\\latexChapters\\bid_sim_covs_sd.tex"), include.rownames = TRUE, 
	  hline.after = hlineAfter,size = "\\setlength{\\tabcolsep}{0.1cm}", 
      table.placement = getOption("xtable.table.placement", paste(table_placement)))


covDistMatrix 		<- lapply(allResults, function(n) n$covDistMatrixLat)
covDistMatrixSD	<- transfBidSimObj(covDistMatrix)

print(xtable(covDistMatrixSD) , type="latex", file=path, include.rownames = TRUE, 
      hline.after = hlineAfter,size = "\\setlength{\\tabcolsep}{0.1cm}", 
      table.placement = getOption("xtable.table.placement", paste(table_placement)))


# Romano Wolf
romWolf				<- lapply(allResults, function(n) n$romWolfDec) # n <- romWolf[[1]]; dgp <- n[[1]]
romWolfMatrix		<- t(sapply(romWolf, function(n) { 
								res <- 	lapply(n, function(dgp){
											dgpTab	<- round(table(factor(dgp, levels = c("Null", "(V)", "Other with (V)", "Other")))/replications*100,1)
										})
								
								resLong <-  c(do.call("cbind",res)) 
								names(resLong) <- rep(c("Null", "(V)", "Other with (V)", "Other"), times = 4)
								return(resLong)
								}))
romWolfMatrix								
latromWolfMatrix 	<- xtable(romWolfMatrix)
align(latromWolfMatrix) <- c("c|", rep(c("c","c","c|"), times = 4))
print(latromWolfMatrix , type="latex", file=paste0("C:\\Users\\LocalAdmin\\Dropbox\\BID DP\\Manuscript\\latexChapters\\bid_sim_romWolf.tex"), include.rownames = TRUE, 
	 table.placement = getOption("xtable.table.placement", paste(table_placement)))


# Selected Basis:
selBasis				<- lapply(allResults, function(n) n$selBasis) # n <- selBasis[[1]]; dgp <- n[[1]]
selBasisMatrix			<- lapply(selBasis, function(n) { 
								res <- 	t(sapply(n, function(dgp){
											dgpTab	<-round(table(factor(dgp, levels = seq(3,maxAllowedBasis,2)))/replications*100,1)
										}))
								return(res)})
				
cNames 					<- colnames(selBasisMatrix[[1]])
rNames 					<- rownames(selBasisMatrix[[1]])

selBasisResMatrix		<- matrix(NA, ncol = length(cNames), nrow = length(rNames)*length(ns))
colnames(selBasisResMatrix) <- cNames
rownames(selBasisResMatrix) <- paste0(rep(ns, each = length(rNames)), rNames)

for(n in ns){
	for(cname in cNames){
		for(rname in rNames){
			selBasisResMatrix[ paste0(n, rname) , cname] <- selBasisMatrix[[as.character(n)]][rname , cname]
		}
	}
}

round(selBasisResMatrix,1)

selBasisResMatrix 	<- xtable(round(selBasisResMatrix,1))
align(selBasisResMatrix) <- c("c", "c","|c|", rep("c", times = 10))
print(selBasisResMatrix , type="latex", file=paste0("C:\\Users\\LocalAdmin\\Dropbox\\BID DP\\Manuscript\\latexChapters\\bid_sim_selBasis.tex"), 
	include.rownames = TRUE, size = "\\setlength{\\tabcolsep}{0.1cm}",  hline.after = seq(4, 4*length(ns), 4),
	 table.placement = getOption("xtable.table.placement", paste(table_placement)))

# Correlations
usedCorrelations	<- lapply(allResults, function(n) n$corSum) # n <- selBasis[[1]]; dgp <- n[[1]]


usedCorrelationsMatrix		<- matrix(NA, ncol = length(rNames), nrow = length(ns))
colnames(usedCorrelationsMatrix) <- rNames
rownames(usedCorrelationsMatrix) <- ns
for(n in ns){
	for(rname in rNames){
		usedCorrelationsMatrix[ as.character(n) , rname] <- round(usedCorrelations[[as.character(n)]][[rname]]["Mean"],2)
	}
}
usedCorrelationsMatrix 	<- xtable(round(usedCorrelationsMatrix,1))
print(usedCorrelationsMatrix , type="latex", file=paste0("C:\\Users\\LocalAdmin\\Dropbox\\BID DP\\Manuscript\\latexChapters\\bid_sim_selCor.tex"), 
	include.rownames = TRUE, size = "\\setlength{\\tabcolsep}{0.1cm}",
	 table.placement = getOption("xtable.table.placement", paste(table_placement)))


