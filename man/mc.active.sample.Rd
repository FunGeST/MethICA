\name{mc.active.sample}
\alias{mc.active.sample}
\title{mc.active.sample}
\description{Determine the most contributing samples for each component}
\usage{mc.active.sample(MC_object, method=c("absolute", "reference"),bval=NULL , MC_contrib_CpG=NULL, number=round(nrow(MC_object$Sample_contrib)*0.1), ref=NULL)}
\arguments{
\item{MC_object}{methylation components object returned by mc.extract}
\item{method}{choose "absolute" to select samples with the greatest contribution to the component (in absolute value), or "reference" to select samples with the strongest methylation difference with respect to reference samples over the most contributing CpG sites}
\item{bval}{methylation (beta-value) matrix}
\item{MC_contrib_CpG}{most contributing CpG sites for each component returned by the mc.active.CpG function}
\item{number}{number of samples selected (default: top 10 percent)}
\item{ref}{list of samples to be used as reference}
}
\value{List of the most contributing samples for each component
}
\author{Lea Meunier}
\keyword{methods}
