\name{mc.annot}
\alias{mc.annot}
\title{mc.annot}
\description{Compute associations between methylation components and sample annotations by uni and multivariate analysis}
\usage{mc.annot(MC_object, annot, selcol=colnames(annot), save=FALSE, output.directory, multi_theshold=0.001)}
\arguments{
\item{MC_object}{methylation components object returned by mc.extract}
\item{annot}{sample annotation table}
\item{selcol}{selection of annotations to be included in univariate analysis (by default all columns will be used)}
\item{save}{if TRUE, save results to the output.directory}
\item{output.directory}{path to save output}
\item{multi_theshold}{p-value threshold to include an annotation in the multivariate analysis}
}
\value{returns two matrices : p-values of univariate analysis and p-values of multivariate analysis
}
\author{Lea Meunier}
\keyword{association}
