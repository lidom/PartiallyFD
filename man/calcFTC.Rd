\name{calcFTC}
\alias{calcFTC}
\title{
FTC Mean and Covariance Estimator of the Paper "Partially Observed Functional Data: The Case
of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
\description{
bla
}
\usage{
calcFTC(fds, comp_dom, basis_seq, base = "fourier", maxBasisLength = 25, 
		basisChoice = "Med", alpha = 0.05, B = 1000,
		f_stat = f_stat2, derFds = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fds}{
	A matrix of observations \eqn{X_i(t_j), i = 1, \ldots, n, j = 1, \ldots, p} where NA are inserted if the function was not observed. 
}
  \item{comp_dom}{
	A numeric vector indicating the complete domain \eqn{[A,B]}.
}
  \item{basis_seq}{
	A numeric vector indicating the basis dimensions to use for Romano Wolf procedure
}
  \item{base}{
	The basis class. "fourier", "monomial", "power" are possible but the latter are not checked.
}
  \item{maxBasisLength}{
	Numeric; the maximum basis allowed.
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
	If \code{FALSE}, the chosen basis will be differentiated and used. If a dataset is supplied, the FTC estimator will be applied with the derivatives supplied. It has the form \eqn{X_i'(t_j), i = 1, \ldots, n, j = 1, \ldots, p} where NA are inserted if the function was not observed. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{krausMean}{
	The pooled mean estimator defined as in Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777–801.
}
	\item{ftcMean}{
	The pooled mean estimator by "Partially Observed Functional Data: The Case of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
	\item{ftcCov}{
	The covariance estimator defined by "Partially Observed Functional Data: The Case of Systematically Missings" by Dominik Liebl and Stefan Rameseder.
}
\item{krausCov}{
	The covariance estimator defined as in Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777–801.
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
	The Romano Wolf Return:
	\item{bla}{
	 Check
	}
}
\item{cor}{
	The empirical correlation of \eqn{d_i} and \eqn{\xi_{i1}}  
}
}
\references{
Kraus (2015): "Components and completion of partially observed functional data". Journal of the Royal Statistical Society 77(4), 777–801.
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

## Application of the ftc Estimator
maxBasisLength 		<- 51		# The Basis Selection Criterion in BIC 
basisSel			<- "Med" 	# The Basis Selection Criterion in BIC 
B					<- 1000 	# Number of Bootstrap Replications
basis_choice		<- "fourier"# The basis where the functions are projected onto
res 				<- calcFTC(fds = bc_fds[["NEG_HT"]], comp_dom = md_dis[[1]],
								 basis_seq = seq(3,maxBasisLength,2), base = basis_choice,
								 maxBasisLength = maxBasisLength, basisChoice = basisSel,
								 alpha = 0.05, B = B, derFds = bc_fds_der[["NEG_HT"]])
str(res, max =1)
## Romano Wolf Decisions
PartiallyFD:::checkFtcHypothesis(res$romWolf$ent)
res$romWolf$pval
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
