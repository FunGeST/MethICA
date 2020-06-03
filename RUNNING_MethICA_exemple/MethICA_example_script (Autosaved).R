### 1. Load data
# Methylation matrix, CpG table & sample annotation table
# function chromatin_feature to add chromatine feature to CpG table 


bval = load.RData('~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/bval.Rdata')
annot = load.RData('~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/LICAFR_annot.RData')
annot= factoall(annot)

output.directory = "~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/output/"
if(!file.exists(output.directory)){
  dir.create(output.directory)
}

# Annotate CpG table with various epigenomic features
#CpG_table = load.RData('~/Google Drive/MethICA/data/CPG_feature_Illumina.RData')
#CpG_feature = chromatin.feature(CpG_table = CpG_table, file_CpG_context = "~/Google Drive/MethICA/data/GSE113405_LIV_ADLT.MethylSeekR.segments.bed", name_col_CpG_context = "CpG_context", file_chrom_state = '~/Google Drive/MethICA/data/cst18_liver.RData', name_col_chrom_state = c("state", "active"), file_CGI = "~/Google Drive/MethICA/data/CGI-based_features_hg19.txt", name_col_CGI ="cgi_feature", file_genes = "~/Google Drive/MethICA/data/Gene-based_features_hg19.txt", name_col_genes = c("gene_name", "gene_feature"), file_replication ="~/Google Drive/MethICA/data/HepG2_replication_domains.RData", name_col_replication = "decile", add_seq_info = TRUE, save = TRUE, output.directory = output.directory)


CpG_feature = load.RData('~/Google Drive/MethICA/RUNNING_MethICA_exemple/LICAFR/CpG_feature.Rdata')

### 2. Perform ICA

# Select most variant CpG sites
NmostVar = 10000 #recommended = 200000
mysd <- apply(bval,1,sd)
sel <- order(mysd,decreasing=T)[1:NmostVar]
bval <- bval[sel,];dim(bval)
CpG_feature <- CpG_feature[rownames(bval),]

# verifier que il y a bien les même CpG dans bval et feature
# verifier qu'il y a bien les même echantillions
all(rownames(annot) == colnames(bval))
all(rownames(bval) == rownames(CpG_feature))


# Perform ICA
# function mc_extract (input: bval matrix | outputs: MC object with two matrices giving the contribution of CpGs and samples to each component)
# Perform n ICA, compute stability and choose the most stable iteration

MC_object <- mc.extract(bval, nb_comp = 20, compute_stability = TRUE, nb_iteration = 20, output.directory = output.directory, save = TRUE)


### 3. Charaterize main methylation changes associated with each component
# Extract the most contributing CpG sites for each MC
MC_contrib_CpG <- mc.active.CpG(MC_object, method = "threshold")


# Extract the most différentially methylated samples for each MC
MC_active_sample = mc.activ.sample(MC_object, method = c("absolute", "reference")[2],bval = bval , MC_contrib_CpG = MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1), ref = grep("N", colnames(bval), value = TRUE))


# Represent methylation changes in most contributing tumors vs. normal samples
# function MC_changes (input: MC object, bval matrix, reference and tumors sample IDs | outputs: graphical representation of methylation changes - equivalent Fig. 1b)

mc.change(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref = grep("N", colnames(bval), value = TRUE), output.directory = output.directory)



### 4. Association of MCw with (epi)genomic characteristics
# Compute and represent enrichment of most contributing CpG sites within specific CGI-based features, chromatin states & methylation context
# one or more functions (mc.cgi, mc.cst, mc.context...) taking as input the MC object and CpG table and giving as outputs enrichment score tables and figures.
# These functions can rely on a same core function computing the enrichment score for any type of annotation.

enrich.CpG.feature(MC_object, MC_contrib_CpG, output.directory = output.directory, CpG_feature = CpG_feature)


### 5. Association of MCs with clinical and molecular features
# Compute the association of each MC vs. annotation table
# function mc.annot (adaptation of the geco.groupsVSannot function) taking as input the MC object and pvalue threshold to keep a feature in multivariate analysis.
# The function will perform univariate analysis, selection of features for multivariate analysis and multivariate analysis.
# Outputs : univariate and multivariate p-value tables.
# We can then give a few examples of downstream representation of significant associations (boxplot, correlation plot etc.) but we will not provide dedicated functions.



Samples_association = mc.annot(MC_object, annot = annot , save = TRUE, output.directory = output.directory, seuil_multi = 0.001)
boxplot(MC_object$Sample_contrib[,"MC13"]~ annot[,"CTNNB1.alt"], col = c("grey30", "grey95"), ylab = "Samples contribution", xlab = "CTNNB1 status", main = "MC13 vs CTNNB1 status")

pvaltab_mol = as.matrix(factoall(Samples_association$pval_uni[,2:ncol(Samples_association$pval_uni)]))
echelle_log = 10^-(seq(0,10, length.out=160)) 
echelle_log[length(echelle_log)] = 0 
echelle_log = as.numeric(echelle_log)
names(echelle_log) = 1:160

tmp_mat = pvaltab_mol
for(i in 1:ncol(pvaltab_mol)){
	for(j in 1:nrow(pvaltab_mol)){
		tmp_mat[j,i] = as.numeric(names(echelle_log[which(echelle_log <(pvaltab_mol[j,i]))])[1])
	}
}
class(tmp_mat)
corrplot::corrplot(tmp_mat, is.corr=FALSE, tl.col = 'black')

pvaltab_mol = as.matrix(factoall(Samples_association$pval_multi))
pvaltab_mol[which(is.na(pvaltab_mol))]=1
echelle_log = 10^-(seq(0,10, length.out=160)) 
echelle_log[length(echelle_log)] = 0 
echelle_log = as.numeric(echelle_log)
names(echelle_log) = 1:160

tmp_mat = pvaltab_mol
for(i in 1:ncol(pvaltab_mol)){
	for(j in 1:nrow(pvaltab_mol)){
		tmp_mat[j,i] = as.numeric(names(echelle_log[which(echelle_log <(pvaltab_mol[j,i]))])[1])
	}
}
class(tmp_mat)
corrplot::corrplot(tmp_mat, is.corr=FALSE, tl.col = 'black')





