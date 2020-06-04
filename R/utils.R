# Rd
# description >> Load a Rdata file into an object
# argument
# item >> file.name >> path of the file to load
# value >> object contained in the Rdata 
# author >> Lea Meunier
# keyword >> Rdata
# end

load.RData = function(file.name) {
  file.real.name = load(file.name, verbose = T)
  return(get(file.real.name))
}

# Rd
# description >> Remove factors in a table, data.frame or matrix
# argument
# item >> d >> table frpm which you want to remove factors
# value >> Same table without factors
# author >> Lea Meunier
# keyword >> table
# end

factoall <- function (d) 
{
    n <- ncol(d)
    for (i in 1:n) {
        if (is.factor(d[, i])) {
            d[, i] <- as.character(d[, i])
            na <- which(is.na(d[, i]))
            num <- suppressWarnings(as.numeric(d[, i]))
            nanum <- which(is.na(num))
            if (length(nanum) == length(na)) {
                d[, i] <- num
            }
        }
    }
    d
}



# Rd
# description >> Annotate positions in a data frame (dfPos) using segments in another data frame (dfSegm)
# argument
# item >> table_Pos >> data frame with positions to annotate
# item >> table_Pos.chrom.col >> chromosome column in table_Pos
# item >> table_Pos.pos.col >> position column in table_Pos
# item >> table_Segm >> data frame with segments to use for annotating dfPos
# item >> table_Segm.chrom.col >> chromosome column in table_Segm
# item >> table_Segm.start.col >> chromosome column in table_Segm
# item >> table_Segm.end.col >> end position column in table_Segm
# item >> cols_to_add >> names of columns in table_Segm that should be used to annotate table_Pos
# item >> names_cols_to_add >> column names to give to the columns added to table_Pos
# value >> table_Pos annotated with table_Segm information
# author >> Lea Meunier
# keyword >> annotation
# end

table.PosXSegm <- function (table_Pos = NULL, table_Pos.chrom.col = "chrom", table_Pos.pos.col = "pos", 
    table_Segm = NULL, table_Segm.chrom.col = "chrom", table_Segm.start.col = "start", 
    table_Segm.end.col = "end", cols_to_add = NULL, names_cols_to_add = NULL) 
{
    for (col in names_cols_to_add){
    	table_Pos[, col] <- NA
    } 
    table_Segm. <- split(table_Segm, table_Segm[, table_Segm.chrom.col])
    table_Pos. <- split(table_Pos, table_Pos[, table_Pos.chrom.col])
    
	for (chr in intersect(names(table_Segm.), names(table_Pos.))) {
		ind <- sapply(table_Pos.[[chr]][, table_Pos.pos.col], 
			function(pos) {
				tmp <- which(table_Segm.[[chr]][, table_Segm.start.col] <= 
				pos & table_Segm.[[chr]][, table_Segm.end.col] >= pos)
				tmp
			})
		for (j in 1:length(cols_to_add)) {
			table_Pos.[[chr]][, names_cols_to_add[j]] <- unlist(lapply(ind, 
				function(z) paste(table_Segm.[[chr]][z, cols_to_add[j]], 
					collapse = ",")))
			table_Pos.[[chr]][which(sapply(ind, length) == 0), names_cols_to_add[j]] <- NA
		}
	}
   
    table_Pos <- unsplit(table_Pos., table_Pos[, table_Pos.chrom.col])
    return(table_Pos)
}


# Rd
# description >> Compute the enrichment of CpG sites within categories of (epi)genomic features
# argument
# item >> CpG_select >> selected CpGs (e.g. most contributing CpGs to a component)
# item >> column >> column in CpG_feature table for which enrichment scores should be computed
# item >> CpG_feature >> annotate CpG tables
# value >> Returns enrichment scores of selected CpGs within each category of the "column" table
# author >> Lea Meunier
# keyword >> test
#` @import stringr
# end

enrich.test <- function (CpG_select, column, CpG_feature) 
{
    total = table(unlist(stringr::str_split(CpG_feature[,column], ",")))
    nom_column = names(total)
    
    selection = table(unlist(stringr::str_split(CpG_feature[CpG_select,column], ",")))[nom_column]
    enrich = (selection/sum(selection, na.rm = TRUE))/(total/sum(total, na.rm = TRUE))
    names(enrich) = nom_column
    return(enrich)
}







