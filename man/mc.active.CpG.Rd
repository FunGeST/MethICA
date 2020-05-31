\name{mc.active.CpG}
\alias{mc.active.CpG}
\title{mc.active.CpG}
\description{Determine the most contributing CpGs for each component}
\usage{mc.active.CpG(MC_object, method=c("threshold", "number"), threshold=quantile(abs(MC_object$CpG_contrib), probs=0.95), number=nrow(MC_object$CpG_contrib)*0.025)}
\arguments{
\item{MC_object}{Object return by mc.extract}
\item{method}{choose between "threshold" selection of CpG above a threshold or number selection of most contributing CpG}
\item{threshold}{threshold value, default select value where 5% of all CpG in all components are higher}
\item{number}{number or most contributingg CpG}
}
\value{return list with active CpG for each component
}
\author{Lea Meunier}
\keyword{methods}
