\name{mc.active.CpG}
\alias{mc.active.CpG}
\title{mc.active.CpG}
\description{Determine the most contributing CpGs for each component}
\usage{mc.active.CpG(MC_object, method=c("threshold", "number"), threshold=quantile(abs(MC_object$CpG_contrib), probs=0.95), number=nrow(MC_object$CpG_contrib)*0.025)}
\arguments{
\item{MC_object}{methylation components object returned by mc.extract}
\item{method}{choose "threshold" to select CpGs with a contribution greater than the "threshold" parameter, or "number" to select a defined number of most contributing CpGs}
\item{threshold}{threshold on contribution values used to select CpGs with the "threshold" method. Set by default to the 95th percentile of CpG contribution values.}
\item{number}{number or most contributing CpGs to select with the "number" method}
}
\value{List of the most contributing CpG sites for each component
}
\author{Lea Meunier}
\keyword{methods}
