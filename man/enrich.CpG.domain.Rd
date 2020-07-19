\name{enrich.CpG.domain}
\alias{enrich.CpG.domain}
\title{enrich.CpG.domain}
\description{Compute and represent enrichment of 48 CpG categories as in Zhou W et al. (Nat Genet 2018)}
\usage{enrich.CpG.domain(CpG_feature, MC_contrib_CpG, MC_active_sample)}
\arguments{
\item{CpG_feature}{table with CpG annotation }
\item{MC_contrib_CpG}{most contributing CpG sites for each component returned by the mc.active.CpG function}
\item{MC_active_sample}{most contributing samples for each component returned by the mc.active.CpG function}
}
\value{return matrice CpG_feature with 48 CpG categories as in Zhou W et al. (Nat Genet 2018)
}
\author{Lea Meunier}
\keyword{methods}
