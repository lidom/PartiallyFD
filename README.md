#### PartiallyFD - A _R_-Package
This readme.md contains an installation and application demo of the code for the Paper "Partially Observed Functional Data: The Case of Systematically Missing" by Dominik Liebl and Stefan Rameseder

##### Installation of required packages
```r
install.packages("devtools")
install.packages("fda")
library("devtools")
```
###### Installation of _PartiallyFD_-Package
```r
install_github("stefanrameseder/PartiallyFD")
library("PartiallyFD")
```
###### Provide Data
```r
data(partObsBidcurves)
attach(partObsBidcurves)
```
Load the combined data for NEG 
```r
combinedNEG <- cbind(bc_fds[["NEG_HT"]], bc_fds[["NEG_NT"]])

matplot(x = md_dis[[1]], y=combinedNEG, ylim = c(0,20000), 
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")
```
Exclude outliers according to Ocker, F., Ehrhart, K.-M., Ott, M. (2015):
_An Economic Analysis of the German Secondary Balancing Power Market, Working Paper (under review)._
```r
firstOutlier          <- which.max(apply(combinedNEG, 2, max, na.rm = TRUE))
secondOutlier         <- which.max(apply(combinedNEG[ ,-firstOutlier], 2, max, na.rm = TRUE))
combinedNEG_woOutlier <- combinedNEG[ ,-c(firstOutlier, secondOutlier+1)]


matplot(x = md_dis[[1]], y=combinedNEG_woOutlier, ylim = c(0,20000), 
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")
```
Load the corresponding data for the derivatives
```r
combinedNEG_der           <- cbind(bc_fds_der[["NEG_HT"]], bc_fds_der[["NEG_NT"]])
combinedNEG_woOutlier_der <- combinedNEG_der[ ,-c(firstOutlier, secondOutlier+1)]
```

##### Application of the ftc Estimator
```r
maxBasisLength    <- 51        # b_max 
basisSel          <- "Med"     # The Basis Selection Criterion over all single BICs
B                 <- 1000      # Number of Bootstrap Replications
basis_choice      <- "fourier" # The basis where the functions are projected onto
derFds            <- combinedNEG_woOutlier_der
res               <- calcFTC(fds = combinedNEG_woOutlier, comp_dom = md_dis[[1]],
                             basis_seq = seq(3,maxBasisLength,2), base = basis_choice,
                             maxBasisLength = maxBasisLength, basisChoice = basisSel,
                             alpha = 0.05, B = B, derFds = combinedNEG_woOutlier_der)
```
Romano Wolf Decision
```r
PartiallyFD:::checkFtcHypothesis(res$romWolf$ent)
```
