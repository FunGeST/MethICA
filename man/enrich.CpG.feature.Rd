\name{enrich.CpG.feature}
\alias{enrich.CpG.feature}
\title{enrich.CpG.feature}
\description{Compute and represent enrichement of the most contributing CpG features to each component within specific (epi)genomic features}
\usage{enrich.CpG.feature(MC_object, MC_contrib_CpG, output.directory, CpG_feature, other_feature_to_test=NA)}
\arguments{
\item{MC_object}{methylation components object returned by mc.extract}
\item{MC_contrib_CpG}{most contributing CpG sites for each component returned by the mc.active.sample function}
\item{output.directory}{path to save output}
\item{other_feature_to_test}{liste of other CpG feature to test  }
}
\author{Lea Meunier}
\keyword{visualisation}
