# Rd
# description >> Annotate a CpG table with various (epi)genomic features
# argument
# item >> CpG_table >> Table giving the genomic position of each CpG. A "FORWARD_SEQUENCE" is required if add_seq_info = TRUE  
# item >> file_CpG_context >> Path of the file containing methylation domains (HMD, PMD, LMR, UMR). Can be a .txt, .bed or .Rdata file with columns c("chr", "start", "end") and columns to add to CpG table
# item >> name_col_CpG_context >> Names of the columns in file_CpG_context that should be used to annotate CpG_table
# item >> file_chrom_state >> Path of the file containing chromatin states. Can be a .txt, .bed or .Rdata file with columns c("chr", "start", "end") and columns to add to CpG table
# item >> name_col_chrom_state >> Names of the columns in file_chrom_state that should be used to annotate CpG_table
# item >> file_CGI >> Path of the file containing CpG island-based features. Can be a .txt, .bed or .Rdata file with columns c("chr", "start", "end") and columns to add to CpG table
# item >> name_col_CGI >> Names of the columns in file_CGI that should be used to annotate CpG_table
# item >> file_genes >> Path of the file containing gene-based features. Can be a .txt, .bed or .Rdata file with columns c("chr", "start", "end") and columns to add to CpG table
# item >> name_col_genes >> Names of the columns in file_genes that should be used to annotate CpG_table
# item >> file_replication >> Path of the file containing replication timing domains. Can be a .txt, .bed or .Rdata file with columns c("chr", "start", "end") and columns to add to CpG table
# item >> name_col_replication >> Names of columns in file_replication that should be used to annotate CpG_table
# item >> add_seq_info >> If add_seq_info == TRUE, add the number of adjacent CpG and the nucleotide context to the CpG_table (cf Zhou et al., Nat Genet 2018)
# value >> returns the CpG table with extra columns indicating the (epi)genomic context of each CpG
# author >> Léa Meunier
# keyword >> CpG annotation
#` @import stringr
# end

chromatin.feature <- function(CpG_table = CpG_table, file_CpG_context = NULL, name_col_CpG_context = NULL, file_chrom_state = NULL, name_col_chrom_state = NULL, file_CGI = NULL, name_col_CGI = NULL, file_genes = NULL, name_col_genes = NULL, file_replication = NULL, name_col_replication = NULL, add_seq_info = TRUE, save = FALSE, output.directory){
	
	CpG_table = data.frame(CpG_table)
	
	## Error if CpG_table does'nt contain minimal column c("TargetID", "MAPINFO", "CHR")
	missing_col = setdiff(c("TargetID", "MAPINFO", "CHR"), colnames(CpG_table))
	if(length(missing_col)>= 1)
		stop(patse0("column", missing_col,"not found in CpG_table"))
	
	# Add several CpG feature
	# Chromatin state
	if(!is.null(file_chrom_state)){
		
		file_name = stringr::str_split(file_chrom_state,"\\.", n = Inf, simplify = TRUE)
		extention = tolower(file_name[ncol(file_name)])
		
		if(extention == "txt" | extention == "bed"){
			chrom_state = read.table(file_chrom_state , sep = "\t", header = T)
		}else if(extention == "rdata"){
			chrom_state = load.RData(file_chrom_state)
		}
		CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = chrom_state, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_chrom_state, names_cols_to_add = c("state", "domain"))
	}
	
	# CpG context
	if(!is.null(file_CpG_context)){
		
		file_name = stringr::str_split(file_CpG_context,"\\.", n = Inf, simplify = TRUE)
		extention = tolower(file_name[ncol(file_name)])
		
		if(extention == "txt" | extention == "bed"){
			CpG_context = read.table(file_CpG_context , sep = "\t", header = T)
			CpG_context$chr = paste("chr", CpG_context$chr, sep = "")
		}else if(extention == "rdata"){
			CpG_context = load.RData(file_CpG_context)
		}
		
		CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = CpG_context, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_CpG_context, names_cols_to_add = "CpG_context")
	}
	
	# CGI feature
	if(!is.null(file_CGI)){
		
		file_name = stringr::str_split(file_CGI,"\\.", n = Inf, simplify = TRUE)
		extention = tolower(file_name[ncol(file_name)])
		
		if(extention == "txt" | extention == "bed"){
			CGI = read.table(file_CGI , sep = "\t", header = T)
		}else if(extention == "rdata"){
			CGI = load.RData(file_CGI)
		}
		
		CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = CGI, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_CGI, names_cols_to_add = "cgi_feature")
	}
	
	# Genes feature
	if(!is.null(file_genes)){
		
		file_name = stringr::str_split(file_genes,"\\.", n = Inf, simplify = TRUE)
		extention = tolower(file_name[ncol(file_name)])
		
		if(extention == "txt" | extention == "bed"){
			Genes = read.table(file_genes , sep = "\t", header = T)
		}else if(extention == "rdata"){
			Genes = load.RData(file_genes)
		}
		
		CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = Genes, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_genes, names_cols_to_add = c("gene_name", "gene_feature"))
	}
	
	# Replication timing
	if(!is.null(file_CGI)){
		
		file_name = stringr::str_split(file_replication,"\\.", n = Inf, simplify = TRUE)
		extention = tolower(file_name[ncol(file_name)])
		
		if(extention == "txt" | extention == "bed"){
			replicatio = read.table(file_replication , sep = "\t", header = T)
		}else if(extention == "rdata"){
			replicatio = load.RData(file_replication)
		}
		
		CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = replicatio, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_replication, names_cols_to_add = "decile")
	}
	
	# Add nucleotide context and number of CpG in the adjacent sequence
	if(add_seq_info == TRUE){
		tmp_table_CpG = data.frame(CpG_table$TargetID , CpG_table$FORWARD_SEQUENCE, str_split(CpG_table$FORWARD_SEQUENCE, "\\[CG\\]", simplify = TRUE))
		colnames(tmp_table_CpG) = c("TargetID", "FORWARD_SEQUENCE", "FORWARD_SEQUENCE_pre", "FORWARD_SEQUENCE_post")
		
		tmp_table_CpG$pre_context = stringr::str_sub(tmp_table_CpG$FORWARD_SEQUENCE_pre, -1, -1)
		tmp_table_CpG$post_context = stringr::str_sub(tmp_table_CpG$FORWARD_SEQUENCE_post, 1, 1)

		CpG_table$context = apply(tmp_table_CpG, 1, function(x){
			if((x[5] == "C" | x[5] == "G")&(x[6] == "C" | x[6] == "G")){
				return("SCGS")
			}else if((x[5] == "C" | x[5] == "G")&(x[6] == "A" | x[6] == "T")){
				return("SCGW")
			}else if((x[5] == "A" | x[5] == "T")&(x[6] == "C" | x[6] == "G")){
				return("SCGW")
			}else{
				return("WCGW")
			}
		})
		
		tmp_table_CpG$FORWARD_SEQUENCE = as.character(tmp_table_CpG$FORWARD_SEQUENCE)
		tmp_table_CpG$FORWARD_SEQUENCE_red = stringr::str_sub(tmp_table_CpG$FORWARD_SEQUENCE, 25, -25)
		tmp_table_CpG$FORWARD_SEQUENCE_red_post = stringr::str_sub(tmp_table_CpG$FORWARD_SEQUENCE_post, 2, -25)

		tmp_table_CpG = data.frame(tmp_table_CpG)
		CpG_table$nb_flanking_CpG = sapply(tmp_table_CpG$FORWARD_SEQUENCE_red, function(x){return(length(gregexpr("CG", x)[[1]])-1)})

	}
	
	if(save == TRUE){
		if(!file.exists(output.directory)){
  			dir.create(output.directory)
		}
		CpG_feature = CpG_table
		save(CpG_feature, file = paste0(output.directory, "CpG_feature.Rdata"))	
		write.table(CpG_feature, file = paste0(output.directory, "CpG_feature.txt"))			
	}
	
	return(CpG_table)
} 
