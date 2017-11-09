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
data(partObsBidcurves)
attach(partObsBidcurves)



## Load the combined data for NEG
combinedNEG <- cbind(bc_fds[["NEG_HT"]], bc_fds[["NEG_NT"]])

matplot(x = md_dis, y=combinedNEG, ylim = c(0,14000), # dim(combinedNEG)
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")

## Exclude outliers according to Ocker, F., Ehrhart, K.-M., Ott, M. (2015):
## "An Economic Analysis of the German Secondary Balancing Power Market?, Working Paper (under review)."
firstOutlier 	<- which.max(apply(combinedNEG, 2, max, na.rm = TRUE))
secondOutlier 	<- which.max(apply(combinedNEG[ ,-firstOutlier], 2, max, na.rm = TRUE))
combinedNEG_woOutlier <- combinedNEG[ ,-c(firstOutlier, secondOutlier+1)]


matplot(x = md_dis, y=combinedNEG_woOutlier, ylim = c(0,20000), # dim(combinedNEG)
        col = "black", type = "l", lwd = 1,
        ylab = "", xaxt = "n", yaxt = "n", lty = "solid")


## Load the corresponding data for the derivatives
combinedNEG_der <- cbind(bc_fds_der[["NEG_HT"]], bc_fds_der[["NEG_NT"]])
combinedNEG_woOutlier_der <- combinedNEG_der[ ,-c(firstOutlier, secondOutlier+1)]
 

## Application of the ftc Estimator
maxBasisLength 		<- 51		# The Basis Selection Criterion in BIC 
basisSel			<- "Med" 	# The Basis Selection Criterion in BIC 
B					<- 1000 	# Number of Bootstrap Replications
basis_choice		<- "fourier"# The basis where the functions are projected onto
derFds 				<- combinedNEG_woOutlier_der
res 				<- calcFTC(fds = combinedNEG_woOutlier, comp_dom = md_dis,
								 basis_seq = seq(3,maxBasisLength,2), base = basis_choice,
								 maxBasisLength = maxBasisLength, basisChoice = basisSel,
								 alpha = 0.05, B = B, derFds = combinedNEG_woOutlier_der)

## Romano Wolf Decisions
PartiallyFD:::checkFtcHypothesis(res$romWolf$ent)
res$romWolf$pval



## #########################################################################
## Both empirical plots as in the paper

maxVals     <- apply(combinedNEG_woOutlier, 2, max, na.rm = TRUE)							 
scl			<- 1.5


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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
