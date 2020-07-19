\name{mc.extract}
\alias{mc.extract}
\title{mc.extract}
\description{Extract Methylation components}
\usage{mc.extract(bval, nb_comp=20, compute_stability=TRUE, nb_iteration=100, alg.typ="parallel", method="C", maxit=1000, fun=c("logcosh","exp"), tol=10^-6, save=FALSE, output.directory)}
\arguments{
\item{bval}{methylation (beta-value) matrix}
\item{nb_comp}{number of components to extract}
\item{compute_stability}{if = TRUE, compute stability of components by performing n iterations and choose the most stable iteration}
\item{nb_iteration}{number of iterations required to assess the stability of components}
\item{alg.typ}{if alg.typ == "parallel" the components are extracted simultaneously (the default). if alg.typ == "deflation" the components are extracted one at a time}
\item{method}{if method == "R" then computations are done exclusively in R (default). The code allows the interested R user to see exactly what the algorithm does. if method == "C" then C code is used to perform most of the computations, which makes the algorithm run faster. During compilation the C code is linked to an optimized BLAS library if present, otherwise stand-alone BLAS routines are compiled.}
\item{maxit}{maximum number of iterations to perform}
\item{fun}{the functional form of the G function used in the approximation to neg-entropy (see ‘details’ of fastICA fonction).}
\item{tol}{a positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged}
\item{save}{if = TRUE, save output in the output.directory}
\item{output.directory}{path to save output}
}
\value{MC_object containing the following elements :
}
\author{Lea Meunier}
\keyword{Independent}
