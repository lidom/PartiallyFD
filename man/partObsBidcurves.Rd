\name{partObsBidcurves}
\alias{partObsBidcurves}
\docType{data}
\alias{partObsBidcurves}
\title{A differentiable functional data set of auction curves}
\description{
  This data set contains weekly functional bidcurves and their derivatives for all four auctions ("POS_HT", "POS_NT", "NEG_HT", "NEG_NT") 
  of Secondary Reserve (SRL) in the German Balance Power Market (www.regelleistung.net).  
}
\usage{data("partObsBidcurves")}
\format{
	A list containing three entries
	 \describe{
    \item{Easy}{
      The beta function of the functional part is a polynomial of order 3:
      \eqn{\beta(t) = -(t-1)^2 + 2} together with two points of impact
      which are located at 0.3 and 0.6 with coefficients -3 and 3.
    }
    \item{Complicated}{
      The beta function of the functional part is a polynomial of order 4:
      \eqn{\beta(t) = -5(t-0.5)^3 - t + 1} together with three points of
      impact are located at 0.3, 0.4 and 0.6 with coefficients -3, 3 and 3.
    }
    \item{NoPoI}{
      The beta function of the functional part is the same as in
      \code{Easy} but the DGP does not include any further points of impact.
    }
    \item{OnlyPoI}{
      Only the points of impact of \code{Easy} are included.
    }
  }
	
	
	, each a list with the names of the auctions "POS_HT", "POS_NT", "NEG_HT", "NEG_NT". 
	Each entry contains for each week between June 27th, 2011, until April, 17th, 2017, (304 weeks in total)
	a differentiable functional object \eqn{X_i(t)}. 
	The functional object is a list of 11 objects where the two functions \code{constrained.curve} and \code{unconstrained.curve}
	can be evaluated using the generic \code{evaluate}. An example can be obtained in the \code{calcFTC} Function.
}
\source{www.regelleistung.net}
\references{
	Liebl, D., and Rameseder, S. (2017): Partially Observed Functional Data: The Case
	of Systematically Missings. Available from the authors on request.
}
\keyword{datasets}
