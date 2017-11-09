### PartiallyFD
This repository contains a R package for the Paper "Partially Observed Functional Data: The Case of Systematically Missing" by Dominik Liebl and Stefan Rameseder. Below, we present a demo for the installation of the package and the application of the estimator. 

#### Installation of required packages
```r
install.packages("devtools")
install.packages("fda")
library("devtools")
```
##### Installation of _PartiallyFD_-Package
```r
install_github("stefanrameseder/PartiallyFD")
library("PartiallyFD")
```
#### Provide Data
```r
data(partObsBidcurves)
attach(partObsBidcurves)
```
###### Load the combined data for NEG 
```r
combinedNEG <- cbind(bc_fds[["NEG_HT"]], bc_fds[["NEG_NT"]])

matplot(x = md_dis[[1]], y=combinedNEG, ylim = c(0,20000), 
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")
```
###### Exclude outliers according to Ocker, F., Ehrhart, K.-M., Ott, M. (2015):_An Economic Analysis of the German Secondary Balancing Power Market, Working Paper (under review)._
```r
firstOutlier          <- which.max(apply(combinedNEG, 2, max, na.rm = TRUE))
secondOutlier         <- which.max(apply(combinedNEG[ ,-firstOutlier], 2, max, na.rm = TRUE))
combinedNEG_woOutlier <- combinedNEG[ ,-c(firstOutlier, secondOutlier+1)]


matplot(x = md_dis[[1]], y=combinedNEG_woOutlier, ylim = c(0,20000), 
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")
```
###### Load the corresponding data for the derivatives
```r
combinedNEG_der           <- cbind(bc_fds_der[["NEG_HT"]], bc_fds_der[["NEG_NT"]])
combinedNEG_woOutlier_der <- combinedNEG_der[ ,-c(firstOutlier, secondOutlier+1)]
```

#### Application of the ftc Estimator
```r
maxBasisLength    <- 51        # b_max 
basisSel          <- "Med"     # The Basis Selection Criterion over all single BICs
B                 <- 1000      # Number of Bootstrap Replications
basis_choice      <- "fourier" # The basis where the functions are projected onto
derFds            <- combinedNEG_woOutlier_der
res               <- calcFTC(fds = combinedNEG_woOutlier, comp_dom = md_dis,
                             basis_seq = seq(3,maxBasisLength,2), base = basis_choice,
                             maxBasisLength = maxBasisLength, basisChoice = basisSel,
                             alpha = 0.05, B = B, derFds = combinedNEG_woOutlier_der)
```
##### Romano Wolf Decision
```r
PartiallyFD:::checkFtcHypothesis(res$romWolf$ent)
```

#### Plot of the Application
```r
maxVals      <- apply(combinedNEG_woOutlier, 2, max, na.rm = TRUE)							 
scl          <- 1.5


layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
par(mar=c(5.1, 5.1, 2.1, 1.1))
matplot(x = 0, y = 0, type = "l", col =PartiallyFD:::addAlpha("black", 0.3), ylim = c(0, 14500), xlim = c(0, 2500), lty = 1, xaxt = "n", yaxt = "n",
        ylab = "Capacity Price in (Euro/MW)/week", xlab = "Electricity Demand in MW", cex.lab=scl+0.45)
for(i in 1:length(maxVals)){ # i = 283
	ind <- which.max(combinedNEG_woOutlier[!is.na(combinedNEG_woOutlier[ ,i]),i])
    if(ind !=1){
		abline(v = md_dis[ind],lwd = 2 , col = "lightblue")
	} 
}
	
abline(v=0, col = "white", lwd = 2)
matplot(x = md_dis, y = combinedNEG_woOutlier, type = "l", col = "black", ylim = c(0, 5000), lty = 1, add = TRUE)

lines(x = md_dis, y = res$ftcMean, col = "blue", lwd = 3)
lines(x = md_dis, y = res$krausMean, col = "darkred", lwd = 3, lty = "dashed")
rect(1750,0,2500,5000,col = rgb(0.5,0.5,0.5,1/4))
axis(side = 1, at = seq(0, 2500, 500),  labels =  paste0(seq(0, 2500, 500), " MW"), cex.axis = scl)	
axis(side = 2, at = seq(0, 14000, 2000),  labels =  paste0(seq(0, 14000, 2000)), cex.axis = scl)	
legend(	"topleft", col = c("black", "darkred", "blue", "lightblue"), inset = 0.01, cex=scl+ 0.2, pt.cex = scl,
        lwd=c(1,3, 3, 3), lty = c("solid", "dashed", "solid","solid"),
        legend = c(expression(paste(X[t])), expression(paste(hat(mu))), expression(paste(hat(mu)[FTC])), expression(paste(d[i]))))
matplot(x = 0, y = 0, type = "l", col = PartiallyFD:::addAlpha("black", 0.3), ylim = c(0, 14500), xlim = c(0, 2500), lty = 1, xaxt = "n", yaxt = "n",
        ylab = "Capacity Price in (Euro/MW)/week", xlab = "Electricity Demand in MW", cex.lab=scl+0.45, add = TRUE)

yLim <- c(0,5000)
xLim <- c(1750, 2500)
par(mar=c(5.1, 3.1, 2.1, 1.5))		
matplot(x = md_dis, y = combinedNEG_woOutlier, type = "l", col = PartiallyFD:::addAlpha("black", 0.3), lty = 1, 
        xaxt = "n", yaxt = "n", ylim = yLim, xlim = xLim,
        ylab = "", xlab = "", cex.lab=scl+0.45)
for(i in 1:length(maxVals)){ # i = 100
    ind <- which.max(combinedNEG_woOutlier[ ,i])
    if(ind !=1){
		abline(v = md_dis[ind],lwd = 2 , col = "lightblue")
	} 
}

lines(x = md_dis, y = res$ftcMean, col = "blue", lwd = 3)
lines(x = md_dis, y = res$krausMean, col = "darkred", lwd = 3, lty = "dashed")
axis(side = 1, at = seq(xLim[1], xLim[2], 250),  labels =  paste0(seq(xLim[1], xLim[2], 250), " MW"), cex.axis = scl)	
axis(side = 2, at = seq(yLim[1], yLim[2], 1000),  labels =  paste0(seq(yLim[1], yLim[2], 1000)), cex.axis = scl)	
matplot(x = md_dis, y = combinedNEG_woOutlier, type = "l", col = PartiallyFD:::addAlpha("black", 0.3), ylim = yLim, xlim = xLim, lty = 1, add = TRUE)
```						 
