\name{calcFTC}
\alias{calcFTC}
\title{
FTC Mean and Covariance Estimator of the Paper "Partially Observed Functional Data: The Case
of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
\description{
The function \code{calcFTC} applies the estimators of equation (1)-(4) in the paper. For details, see article.
}
\usage{
calcFTC(fds, comp_dom, basis_seq, base = "fourier", maxBasisLength = 25, 
		basisChoice = "Med", alpha = 0.05, B = 1000,
		f_stat = f_stat2, derFds = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fds}{
	A matrix of observations \eqn{X_i(t_j), i = 1, \ldots, n, j = 1, \ldots, p} where NA are inserted if function \eqn{X_i(\cdot)} was not observed at point \eqn{t_j}. 
}
  \item{comp_dom}{
	A numeric vector indicating the grid on the complete domain \eqn{[A,B]}.
}
  \item{basis_seq}{
	A numeric vector indicating the basis dimensions to use for Romano Wolf procedure. Typically, this is \eqn{3, 5, \ldots, b_m}.
}
  \item{base}{
	The basis class. "fourier", "monomial", "power" are possible but the latter are not checked.
}
  \item{maxBasisLength}{
	Numeric; the maximum basis \eqn{b_m} allowed.
}
  \item{basisChoice}{
    The selection criterion over all \eqn{\hat b_i}. "Max" or "Med".
}
  \item{alpha}{
	The significance level for Romano Wolf.
}
  \item{B}{
	The number of bootstrap replications for Romano Wolf.
}
  \item{f_stat}{
	The functional test statistic which are obtained from the regression of \eqn{d} onto the coefficients \eqn{xi_j}
}
  \item{derFds}{
	If \code{FALSE}, the chosen basis will be differentiated and used. If a dataset is supplied, the FTC estimator will be applied with the derivatives supplied. It has the form \eqn{X_i'(t_j), i = 1, \ldots, n, j = 1, \ldots, p} where \code{NA}s are inserted if the function was not observed. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{krausMean}{
	The pooled mean estimator defined as in Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777?801.
}
	\item{ftcMean}{
	The pooled mean estimator by "Partially Observed Functional Data: The Case of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
	\item{ftcCov}{
	The covariance estimator defined by "Partially Observed Functional Data: The Case of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
\item{krausCov}{
	The covariance estimator defined as in Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777?801.
}
\item{bicMat}{
	A matrix with the different BICs over all bases dimensions.
}
\item{coefs}{
	The coeffiencts of the representation of the \eqn{X_i(t)} in the basis.
}
\item{selBasis}{
	The selected dimension of the basis.
}
\item{romWolf}{
	The Romano Wolf Return
}
\item{cor}{
	The empirical correlation of \eqn{d_i} and \eqn{\xi_{i1}}  
}
}
\references{
Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777?801.
}
\author{
Stefan Rameseder
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(log_partObsBidcurves)
attach(log_partObsBidcurves)

## Load the combined data for NEG
combinedNEG <- cbind(log_bc_fds[["NEG_HT"]], log_bc_fds[["NEG_NT"]])

matplot(x = md_dis, y=combinedNEG, ylim = c(0,10), # dim(combinedNEG)
        col = PartiallyFD:::addAlpha("black", 0.2), type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")

## Exclude outliers according to Ocker, F., Ehrhart, K.-M., Ott, M. (2015):
## "An Economic Analysis of the German Secondary Balancing Power Market, Working Paper (under review)."
firstOutlier 	<- which.max(apply(combinedNEG, 2, max, na.rm = TRUE))
secondOutlier 	<- which.max(apply(combinedNEG[ ,-firstOutlier], 2, max, na.rm = TRUE))
combinedNEG_woOutlier <- combinedNEG[ ,-c(firstOutlier, secondOutlier+1)]

matplot(x = md_dis, y=combinedNEG_woOutlier, ylim = c(0,10), # dim(combinedNEG)
        col = PartiallyFD:::addAlpha("black", 0.2), type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")

## Observe point mass at 0 and mixture of two normals by calculating 
## Principal Components and their Scores
## Find components on "stable" domain [400MW, 1200MW)
d_max       <- which.min(md_dis < 1200)
d_min       <- which.max(md_dis >= 400)

## Calculate Eigenvectors and Scores
X_mat       <- combinedNEG_woOutlier[d_min:d_max, ]
X_cent_mat  <- X_mat - rowMeans(X_mat)
Cov_mat     <- X_cent_mat \%*\% t(X_cent_mat) / ncol(X_cent_mat)
eigen_obj   <- eigen(Cov_mat)

## First five cumulative Variances
c(cumsum(eigen_obj$values)/sum(eigen_obj$values))[1:5]

## Standardizing to L2-norm == 1 (assuming [a,b]=[0,1])
eigenvec_1  <- eigen_obj$vectors[,1] * sqrt(length(d_min:d_max))
eigenvec_1  <- eigenvec_1 * sign(sum(eigenvec_1))
PC_scores   <- c(eigenvec_1 \%*\% X_cent_mat)/length(d_min:d_max)


## Remove point-mass on minimal PC-scores (correspond to zero-functions) 
quantile(PC_scores, seq(0, 0.1, 0.005))
thr         <- -5.16
scoreIndicesProbMass <- PC_scores <= thr; 
PC_scores_red <- PC_scores[PC_scores > thr]

##
mclust.obj <- densityMclust(data = PC_scores_red, G=2) 
clust_vec  <- mclust.obj$classification

## Cluster Plot
par(mfrow=c(1,1), mar=c(4.5,4,2.5,1)+0.1, family = "sans")
plotDensityMclust1(mclust.obj, xlab="First FPC-Scores (99.8\% Explained Variance)", main="")
mtext(text = "Normal Mixture Cluster Result", side = 3, line = 1.25, cex=1.2)
hist(PC_scores_red, add=TRUE, freq = FALSE, breaks = 12)
points(x =              PC_scores_red[clust_vec==2], 
       y = rep(0,length(PC_scores_red[clust_vec==2])), col = PartiallyFD:::addAlpha("red",1), bg=PartiallyFD:::addAlpha("red",1),
       pch=21, cex=1)
points(x =              PC_scores_red[clust_vec==1], 
       y = rep(0,length(PC_scores_red[clust_vec==1])), col = PartiallyFD:::addAlpha("blue",1), bg=PartiallyFD:::addAlpha("blue",1), 
       pch=22, cex=1)
points(x = PC_scores[PC_scores < thr], 
       y = seq(0,0.1,len=length(PC_scores[PC_scores < thr])), col=PartiallyFD:::addAlpha("darkorange",1), bg=PartiallyFD:::addAlpha("darkorange",1),
       pch=23, cex=1)
##
legend("topleft", legend = c("High-Price Cluster", "Low-Price Cluster", "Zero-Functions"), 
       pt.cex=1,
       pch=c(21,22,23), 
       pt.bg = c(PartiallyFD:::addAlpha("red",1), PartiallyFD:::addAlpha("blue",1), PartiallyFD:::addAlpha("darkorange",1)), 
       col=c(PartiallyFD:::addAlpha("red",1), PartiallyFD:::addAlpha("blue",1), PartiallyFD:::addAlpha("darkorange",1)), bty="n")
dev.off()

## Testing normality using the KS-Test:
ks.test(PC_scores_red[clust_vec==2], "pnorm", mean(PC_scores_red[clust_vec==2]), sd(PC_scores_red[clust_vec==2]))

## Vector of indices of the 'Low-Price Cluster': these will be excluded
scoreIndicesClust <- clust_vec==1
paste0("Point mass at Zero: ",sum(scoreIndicesProbMass) ," and Low-Price Cluster:", sum(scoreIndicesClust))

## Reduce Functions and Derivatives by Low-Price Cluster and Point Mass Curves 
reducedDomSample <- combinedNEG_woOutlier[ , !scoreIndicesProbMass]
combinedNEG_woOutlier <- reducedDomSample[ , !scoreIndicesClust]

matplot(combinedNEG_woOutlier, type = "l", col = PartiallyFD:::addAlpha("black", 0.3))

## Reduce Derivates by Outliers, Low-Price Cluster, and Point Mass Curves
combinedNEG_der <- cbind(log_bc_fds_der[["NEG_HT"]], log_bc_fds_der[["NEG_NT"]])
combinedNEG_woOutlier_der <- combinedNEG_der[ , -c(firstOutlier, secondOutlier+1)]
reducedDomSample <- combinedNEG_woOutlier_der[ , !scoreIndicesProbMass]
combinedNEG_woOutlier_der <- reducedDomSample[ , !scoreIndicesClust]

## Application of the ftc Estimator
maxBasisLength 		<- 51		# The Basis Selection Criterion in BIC 
basisSel			<- "Med" 	# The Basis Selection Criterion in BIC 
B					<- 1000 	# Number of Bootstrap Replications
basis_choice		<- "fourier"# The basis where the functions are projected onto
res 				<- calcFTC(fds = combinedNEG_woOutlier, comp_dom = md_dis,
                   basis_seq = seq(3,maxBasisLength,2), base = basis_choice,
                   maxBasisLength = maxBasisLength, basisChoice = basisSel,
                   alpha = 0.05, B = B, derFds = combinedNEG_woOutlier_der)

## Romano Wolf Decisions
PartiallyFD:::checkFtcHypothesis(res$romWolf$ent)

## Calculate Confidence Intervals
alpha           <- 0.05
nPerP           <- apply(combinedNEG_woOutlier, 1, function(x) sum(!is.na(x))) 
ftcSD           <- sqrt(diag(res$ftcCov))/sqrt(nPerP)
critValue       <- qnorm(1-(alpha/2))
CISummand       <- critValue * ftcSD
ftcCI_plus      <- res$ftcMean + CISummand
ftcCI_minus     <- res$ftcMean - CISummand

# Kraus Confidence Intervalls
krausSD         <- sqrt(diag(res$krausCov))/sqrt(nPerP)
CISummand       <- critValue * krausSD
krausCI_plus    <- res$krausMean + CISummand
krausCI_minus   <- res$krausMean - CISummand

## Both empirical plots as in the paper
scl.axs             <- 1.9
p                   <- length(res$krausMean)
p_seq               <- seq(1,p,8)
p.cex               <- 1.2
maxVals     		<- apply(combinedNEG_woOutlier, 2, max, na.rm = TRUE)							 


layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
par(mar=c(5.1, 5.1, 2.1, 1.1))
matplot(x = 0, y = 0, type = "l", col = PartiallyFD:::addAlpha("black", 0.2), ylim = c(0, 10), xlim = c(0, 2500), lty = 1, xaxt = "n", yaxt = "n",
        ylab = "Log Price [Log(Euro/MW)/week]", xlab = "Electricity Demand [MW]", cex.lab=scl.axs+0.45)

for(i in 1:length(maxVals)){ # i = 283
	ind <- which.max(combinedNEG_woOutlier[!is.na(combinedNEG_woOutlier[ ,i]),i])
    if(ind !=1){
		abline(v = md_dis[ind],lwd = 2 , col = "lightblue")
	} 
}
matplot(x = md_dis, y = combinedNEG_woOutlier, type = "l", col = PartiallyFD:::addAlpha("black", 0.2), ylim = c(0, 10), lty = 1, add = TRUE)

points(x = md_dis[p_seq], y = res$ftcMean[p_seq], col = "blue", pch = 16, cex = p.cex)
lines(x = md_dis, y = res$ftcMean, col = "blue", lwd = 2)
points(x = md_dis[p_seq], y = res$krausMean[p_seq], col = "darkred", pch = 18, cex = p.cex)
lines(x = md_dis, y = res$krausMean, col = "darkred", lwd = 2)

polygon(x = c(md_dis, rev(md_dis)),
        y = c(ftcCI_plus, rev(ftcCI_minus)),
        col = PartiallyFD:::addAlpha("blue", 0.2), border = "blue", lwd = 1)
polygon(x = c(md_dis, rev(md_dis)),
        y = c(krausCI_plus, rev(krausCI_minus)),
        col = PartiallyFD:::addAlpha("darkred", 0.2), border = "darkred", lwd = 1)

yLim <- c(4,10)
xLim <- c(1750, 2500)

rect(xLim[1],yLim[1],xLim[2],yLim[2],col = rgb(0.5,0.5,0.5,1/4))
axis(side = 1, at = seq(0, 2500, 500),  labels =  paste0(seq(0, 2500, 500), " MW"), cex.axis = scl.axs)	
axis(side = 2, at = seq(0, 10, 2),  labels =  paste0(seq(0, 10, 2)), cex.axis = scl.axs)	
legend(	"topleft", col = c("black", "darkred", "blue", "lightblue"), inset = 0.01, cex=scl.axs+ 0.4, pt.cex = scl.axs,
        lwd=c(1,2, 2, 3), lty = c("solid", "solid", "solid","solid"),
        pch=c(NA,18,16, NA),  
        legend = c(expression(paste(X[t])), expression(paste(hat(mu))), expression(paste(hat(mu)[FTC])), expression(paste(d[i]))))
matplot(x = 0, y = 0, type = "l", col = PartiallyFD:::addAlpha("black", 0.2), ylim = c(0, 10), xlim = c(0, 2500), lty = 1, xaxt = "n", yaxt = "n",
        ylab = "Log Price in (Log Euro/MW)/week", xlab = "Electricity Demand in MW", cex.lab=scl.axs+0.45, add = TRUE)


par(mar=c(5.1, 3.1, 2.1, 2.1))		
matplot(x = md_dis, y = combinedNEG, type = "l", col = PartiallyFD:::addAlpha("black", 0.2), lty = 1, 
        xaxt = "n", yaxt = "n", ylim = yLim, xlim = xLim,
        ylab = "", xlab = "", cex.lab=scl.axs+0.45)
for(i in 1:length(maxVals)){ # i = 283
	ind <- which.max(combinedNEG_woOutlier[!is.na(combinedNEG_woOutlier[ ,i]),i])
    if(ind !=1){
		abline(v = md_dis[ind],lwd = 2 , col = "lightblue")
	} 
}

points(x = md_dis[p_seq], y = res$ftcMean[p_seq], col = "blue", pch = 16, cex = p.cex)
lines(x = md_dis, y = res$ftcMean, col = "blue", lwd = 2)
points(x = md_dis[p_seq], y = res$krausMean[p_seq], col = "darkred", pch = 18, cex = p.cex)
lines(x = md_dis, y = res$krausMean, col = "darkred", lwd = 2)

polygon(x = c(md_dis, rev(md_dis)),
        y = c(ftcCI_plus, rev(ftcCI_minus)),
        col = PartiallyFD:::addAlpha("blue", 0.2), border = "blue", lwd = 1)
polygon(x = c(md_dis, rev(md_dis)),
        y = c(krausCI_plus, rev(krausCI_minus)),
        col = PartiallyFD:::addAlpha("darkred", 0.2), border = "darkred", lwd = 1)
axis(side = 1, at = seq(xLim[1], xLim[2], 250),  labels =  paste0(seq(xLim[1], xLim[2], 250), " MW"), cex.axis = scl.axs)	
axis(side = 2, at = seq(yLim[1], yLim[2], 2),  labels =  paste0(seq(yLim[1], yLim[2], 2)), cex.axis = scl.axs)	
matplot(x = md_dis, y = combinedNEG, type = "l", col = PartiallyFD:::addAlpha("black", 0.3), ylim = yLim, xlim = xLim, lty = 1, add = TRUE)
dev.off()
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
