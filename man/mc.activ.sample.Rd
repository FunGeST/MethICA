\name{mc.activ.sample}
\alias{mc.activ.sample}
\title{mc.activ.sample}
\description{Determine the most active sample in each component}
\usage{mc.activ.sample(MC_object, method=c("absolute", "reference"),bval=NULL , MC_contrib_CpG=NULL, number=round(nrow(MC_object$Sample_contrib)*0.1), ref=NULL)}
\arguments{
\item{MC_object}{Object return by mc.extract}
\item{method}{choose between "absolute" selection of most contributing Samples in the component (in absolute value) and "reference", the samples with bval of active CpG farthest for bval in ref samples}
\item{bval}{bval methylation data}
\item{MC_contrib_CpG}{most contributing CpG of each comp compute with mc.active.CpG function}
\item{number}{nomber of samples selected (default : 10 percent)}
\item{ref}{liste of samples use like reference}
}
\value{return liste of active sample in each component
}
\author{Lea Meunier}
\keyword{methods}
