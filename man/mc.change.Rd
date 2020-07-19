\name{mc.change}
\alias{mc.change}
\title{mc.change}
\description{Generate graphs showing methylation changes between reference samples and the most contributing samples to each component}
\usage{mc.change(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref, output.directory)}
\arguments{
\item{MC_object}{methylation components object returned by mc.extract}
\item{MC_active_sample}{most contributing samples for each component returned by the mc.active.CpG function}
\item{MC_contrib_CpG}{most contributing CpG sites for each component returned by the mc.active.sample function}
\item{bval}{methylation (beta-value) matrix}
\item{ref}{list of samples to be used as reference}
\item{output.directory}{path to save output}
}
\author{Lea Meunier}
\keyword{visualisation}
