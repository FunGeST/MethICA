library(MethICA)




### 1. Load data
# Methylation matrix, CpG table & sample annotation table.
test.directory <- "~/Tools/MethICA/RUNNING_MethICA_example/LICAFR"

bval = load.RData('~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/bval.Rdata')

annot = load.RData(file.path(test.directory,'LICAFR_annot.RData'))
annot = factoall(annot)

CpG_feature = load.RData('~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/CpG_feature.Rdata')
# The chromatin_feature function may be used to annotate a CpG table with chromatin features.
#CpG_table = load.RData('~/Google Drive/MethICA/data/CPG_feature_Illumina.RData')
#CpG_feature = chromatin.feature(CpG_table = CpG_table, file_CpG_context = "~/Google Drive/MethICA/data/GSE113405_LIV_ADLT.MethylSeekR.segments.bed", name_col_CpG_context = "CpG_context", file_chrom_state = '~/Google Drive/MethICA/data/cst18_liver.RData', name_col_chrom_state = c("state", "active"), file_CGI = "~/Google Drive/MethICA/data/CGI-based_features_hg19.txt", name_col_CGI ="cgi_feature", file_genes = "~/Google Drive/MethICA/data/Gene-based_features_hg19.txt", name_col_genes = c("gene_name", "gene_feature"), file_replication ="~/Google Drive/MethICA/data/HepG2_replication_domains.RData", name_col_replication = "decile", add_seq_info = TRUE, save = TRUE, output.directory = output.directory)

# Check that CpGs and samples have the same order in bval and annotation tables
all(rownames(annot) == colnames(bval))
all(rownames(bval) == rownames(CpG_feature))

output.directory = file.path(test.directory,"output")
if(!file.exists(output.directory)){ dir.create(output.directory) }




### 2. Perform ICA
# Select most variant CpG sites
NmostVar = 10000 #recommended = 200000
mysd <- apply(bval,1,sd)
sel <- order(mysd,decreasing=T)[1:NmostVar]
bval <- bval[sel,];dim(bval)
CpG_feature <- CpG_feature[rownames(bval),]

# Run ICA
# function mc_extract (input: bval matrix | outputs: MC object with two matrices giving the contribution of CpGs and samples to each component)
# Performs n iterations of ICA, computes stability and chooses the most stable iteration
MC_object <- mc.extract(bval, nb_comp = 20, compute_stability = TRUE, nb_iteration = 20, output.directory = output.directory, save = TRUE)




### 3. Charaterize the main methylation changes associated with each methylation component (MC)
# Extract the most contributing CpG sites for each MC
MC_contrib_CpG <- mc.active.CpG(MC_object, method = "threshold")

# Extract the most contributing samples for each MC
MC_active_sample = mc.active.sample(MC_object, method = c("absolute", "reference")[2],bval = bval , MC_contrib_CpG = MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1), ref = grep("N", colnames(bval), value = TRUE))

# Represent methylation changes in the most contributing samples vs. reference samples
mc.change(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref = grep("N", colnames(bval), value = TRUE), output.directory = output.directory)




### 4. Association of MCs with (epi)genomic characteristics
# Compute and represent enrichment of most contributing CpG sites within specific CGI-based features, chromatin states & methylation contexts.
enrich.CpG.feature(MC_object, MC_contrib_CpG, output.directory = output.directory, CpG_feature = CpG_feature)




### 5. Association of MCs with sample annotations
sample.assoc = mc.annot(MC_object, annot = annot , save = TRUE, output.directory = output.directory, seuil_multi = 0.001)

#Examples of representations for sample associations
#boxplot
boxplot(MC_object$Sample_contrib[,"MC13"]~ annot[,"CTNNB1.alt"], col = c("grey30", "grey95"), ylab = "Sample contribution", xlab = "CTNNB1 status", main = "MC13 vs CTNNB1 status")


#corrplot for univariate analysis
pvaltab_mol = as.matrix(factoall(sample.assoc$pval_uni[,2:ncol(sample.assoc$pval_uni)]))
#need to convert p-value in entier number for the corrplot function
log_scale = 10^-(seq(0,10, length.out=160)) 
log_scale[length(log_scale)] = 0 
log_scale = as.numeric(log_scale)
#scale for the corrplot function
names(log_scale) = 1:160

#assign each p-value to the scale
tmp_mat = pvaltab_mol
for(i in 1:ncol(pvaltab_mol)){
	for(j in 1:nrow(pvaltab_mol)){
		tmp_mat[j,i] = as.numeric(names(log_scale[which(log_scale <(pvaltab_mol[j,i]))])[1])
	}
}
class(tmp_mat)
#use the function 
corrplot::corrplot(tmp_mat, is.corr=FALSE, tl.col = 'black')

#same for multivariate analysis
pvaltab_mol = as.matrix(factoall(sample.assoc$pval_multi))
#replace NA by 1
pvaltab_mol[which(is.na(pvaltab_mol))]=1
log_scale = 10^-(seq(0,10, length.out=160)) 
log_scale[length(log_scale)] = 0 
log_scale = as.numeric(log_scale)
names(log_scale) = 1:160

tmp_mat = pvaltab_mol
for(i in 1:ncol(pvaltab_mol)){
	for(j in 1:nrow(pvaltab_mol)){
		tmp_mat[j,i] = as.numeric(names(log_scale[which(log_scale <(pvaltab_mol[j,i]))])[1])
	}
}
class(tmp_mat)
corrplot::corrplot(tmp_mat, is.corr=FALSE, tl.col = 'black')





