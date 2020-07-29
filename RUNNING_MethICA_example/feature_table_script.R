output.directory = "~/MethICA_feature_annotation"

# Please use this link to downlad test data: https://drive.google.com/drive/folders/1BTQOhvI_qQou1CD94N_TCV_TEbcBC671?usp=sharing
# These files contain various types of (epi)genomic annotations for liver tissue
# Indicate the path of the downloaded folder below.
data.directory <- "~/Downloads/MethICAdata/" # Indicate here the path to test data directory
load(file.path(data.directory,'hg19_genome_feature.Rdata'),verbose = T)
load(file.path(data.directory,'liver_specific_data.Rdata'),verbose = T)
	
CpG_table = data.frame(CpG_table)
	
## CpG_table must contain minimal column c("TargetID", "MAPINFO", "CHR")

# typical script to add several CpG feature
# - import segment table with feature of interest with minimal column "chr", "start", "end"
# - use table.PosXSegm function with the names of column to add, and the names give in the final CpG_feature table

# These annotations take a while. Please be patient!

# Chromatin state
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = cst18, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = c("state", "active"), names_cols_to_add = c("state", "domain"))

# CpG context
GSE113405_LIV_ADLT.MethylSeekR.segments$chr <- paste0("chr",GSE113405_LIV_ADLT.MethylSeekR.segments$chr)
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = GSE113405_LIV_ADLT.MethylSeekR.segments, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = "CpG_context", names_cols_to_add = "CpG_context")

# CGI feature
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = CGIbased_features_hg19, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = "cgi_feature", names_cols_to_add = "cgi_feature")

# Genes feature
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = Genebased_features_hg19, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = c("gene_name", "gene_feature"), names_cols_to_add = c("gene_name", "gene_feature"))
	
# Replication timing
CpG_table = table.PosXSegm(table_Pos = CpG_table, table_Pos.chrom.col = "CHR", table_Pos.pos.col = "MAPINFO", 
    		table_Segm = repseq, table_Segm.chrom.col = "chr", table_Segm.start.col = "start", 
    		table_Segm.end.col = "end", cols_to_add = "decile", names_cols_to_add = "decile")
	
# Add nucleotide context and number of CpG in the adjacent sequence
tmp_table_CpG = data.frame(CpG_table$TargetID , CpG_table$FORWARD_SEQUENCE, stringr::str_split(CpG_table$FORWARD_SEQUENCE, "\\[CG\\]", simplify = TRUE))
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

