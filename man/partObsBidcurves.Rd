\name{partObsBidcurves}
\alias{partObsBidcurves}
\docType{data}
\alias{partObsBidcurves}
\title{A differentiable functional data set of auction curves}
\description{
  This data set contains weekly bidcurves \eqn{X_i(t_j)} and their derivatives \eqn{X_i'(t_j)} for all four auctions ("POS_HT", "POS_NT", "NEG_HT", "NEG_NT") for \eqn{i=1, \ldots, 304} and \eqn{j=1, \ldots, 501} of Secondary Reserve (SRL) in the German Balance Power Market (www.regelleistung.net).  The \eqn{i} correspond to dates June 27th, 2011, until April, 17th, 2017, (304 weeks in total).
}
\usage{data("partObsBidcurves")}
\format{
	A list containing three entries
	 \describe{
    \item{bc_fds}{
      A list of observations \eqn{X_i(t_j)} with four entries for each auction "POS_HT", ..., "NEG_NT" and NAs if the function was not observed. 
    }
    \item{bc_fds_der}{
      A list of observations \eqn{X_i'(t_j)} with four entries for each auction "POS_HT", ..., "NEG_NT" and NAs if the derivative was not observed. 
    }
    \item{md_dis}{
      The grid \eqn{t_1, \ldots, t_p}, where the functions are observed.
    }
  }
}
\source{www.regelleistung.net}
\references{
	Liebl, D., and Rameseder, S. (2017): Partially Observed Functional Data: The Case
	of Systematically Missings. Available from the authors on request.
}
\keyword{datasets}
