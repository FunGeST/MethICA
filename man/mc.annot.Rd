\name{mc.annot}
\alias{mc.annot}
\title{mc.annot}
\description{compute uni and multivariate analysis of component and annotations}
\usage{mc.annot(MC_object, annot, selcol=colnames(annot), save=FALSE, output.directory, selcol_Multi=selcol, seuil_multi=0.001)}
\arguments{
\item{MC_object}{Object return by mc.extract}
\item{annot}{annotation table}
\item{selcol}{selection of annotation for univariate analysis}
\item{save}{if TRUE, save results file in the output.directory}
\item{output.directory}{path to save output}
\item{selcol_Multi}{selection of annotation for multivariate analysis}
}
\value{return 2 matrix : pvalue of univariate analysis and pvalue of multivariate analysis
}
\author{Lea Meunier}
\keyword{association}
