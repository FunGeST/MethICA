\name{chromatin.feature}
\alias{chromatin.feature}
\title{chromatin.feature}
\description{required}
\usage{chromatin.feature(CpG_table=CpG_table, file_CpG_context=NULL, name_col_CpG_context=NULL, file_chrom_state=NULL, name_col_chrom_state=NULL, file_CGI=NULL, name_col_CGI=NULL, file_genes=NULL, name_col_genes=NULL, file_replication=NULL, name_col_replication=NULL, add_seq_info=TRUE, save=FALSE, output.directory)}
\arguments{
\item{CpG_table}{with c("TargetID", "MAPINFO", "CHR") and "FORWARD_SEQUENCE" if add_seq_info = TRUE  }
\item{file_CpG_context}{Path of the file containing CpG context (HMD, PMD, LMR, UMR) file in .txt, .bed or Rdata with CpG feature to add to the table, with column c("chr", "start", "end") and the column to add}
\item{name_col_CpG_context}{Names of columns in file_CpG_context that should be used to annotate CpG_table}
\item{file_chrom_state}{Path of the file containing chromatin states file in .txt, .bed or Rdata with CpG feature to add to the table, with column c("chr", "start", "end") and the column to add}
\item{name_col_chrom_state}{Names of columns in file_chrom_state that should be used to annotate CpG_table}
\item{file_CGI}{Path of the file containing. file in .txt, .bed or Rdata with CpG feature to add to the table, with column c("chr", "start", "end") and the column to add}
\item{name_col_CGI}{Names of columns in file_CGI that should be used to annotate CpG_table}
\item{file_genes}{Path of the file containing. file in .txt, .bed or Rdata with CpG feature to add to the table, with column c("chr", "start", "end") and the column to add}
\item{name_col_genes}{Names of columns in file_genes that should be used to annotate CpG_table}
\item{file_replication}{Path of the file containing. file in .txt, .bed or Rdata with CpG feature to add to the table, with column c("chr", "start", "end") and the column to add}
\item{name_col_replication}{Names of columns in file_replication that should be used to annotate CpG_table}
\item{add_seq_info}{If = TRUE, add the number of adjacent CpG and the nucleotide context to the CpG_table}
}
\value{return CpG_table annoted with Chromatine feature
}
\author{Léa Meunier}
\keyword{CpG}
