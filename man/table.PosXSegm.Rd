\name{table.PosXSegm}
\alias{table.PosXSegm}
\title{table.PosXSegm}
\description{Function to annotate positions in a data frame (dfPos) using segments in another data frame (dfSegm)}
\usage{table.PosXSegm(table_Pos=NULL, table_Pos.chrom.col="chrom", table_Pos.pos.col="pos", table_Segm=NULL, table_Segm.chrom.col="chrom", table_Segm.start.col="start", table_Segm.end.col="end", cols_to_add=NULL, names_cols_to_add=NULL) }
\arguments{
\item{table_Pos}{Data frame with positions to annotate}
\item{table_Pos.chrom.col}{Chromosome column in table_Pos}
\item{table_Pos.pos.col}{Position column in table_Pos}
\item{table_Segm}{required}
\item{table_Segm.chrom.col}{Data frame with segments to use for annotating}
\item{table_Segm.start.col}{Chromosome column in table_Segm}
\item{table_Segm.end.col}{Position column in table_Segm}
\item{cols_to_add}{Names of columns in table_Segm that should be used to annotate table_Pos}
\item{names_cols_to_add}{Column names to give in the columns added to table_Pos}
}
\value{table_Pos annotated with table_Segm information
}
\author{Lea Meunier}
\keyword{annotation}
