output.directory = ./output
	
load(file.path(data.directory,'/data/hg19_genome_feature.Rdata'),verbose = T)
load(file.path(data.directory,'/data/liver_specific_data.Rdata'),verbose = T)
	
CpG_table = data.frame(CpG_table)
	
## CpG_table must contain minimal column c("TargetID", "MAPINFO", "CHR")


# typical script to add several CpG feature
# - import segment table with feature of interest with minimal column "chr", "start", "end"
# - use table.PosXSegm function with the names of column to add, and the names give in the final CpG_feature table


# Chromatin state
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = chrom_state, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_chrom_state, names_cols_to_add = c("state", "domain"))

# CpG context
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = CpG_context, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_CpG_context, names_cols_to_add = "CpG_context")

	
# CGI feature
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = CGI, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_CGI, names_cols_to_add = "cgi_feature")

# Genes feature
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = Genes, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_genes, names_cols_to_add = c("gene_name", "gene_feature"))
	
# Replication timing
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = replicatio, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = name_col_replication, names_cols_to_add = "decile")

	
# Add nucleotide context and number of CpG in the adjacent sequence
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


CpG_feature = CpG_table
save(CpG_feature, file = file.path(output.directory, "CpG_feature.Rdata"))	
write.table(CpG_feature, file = file.path(output.directory, "CpG_feature.txt"))			

