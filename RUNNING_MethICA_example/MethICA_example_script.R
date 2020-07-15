library(MethICA)
library(MethICAdata)
library(corrplot)

data.directory <- file.path(.libPaths(), MethICAdata)


### 1. Load data
# Methylation matrix, CpG table & sample annotation table.
output.directory = "~/Results/"
if(!file.exists(output.directory)){ dir.create(output.directory) }

data.directory <- file.path(.libPaths(), MethICAdata)
load(file.path(data.directory,'/data/LICAFR_methylation.Rdata'),verbose = T)

# To create the file CpG_feature use in the analysis, see the feature_table_script.R

### 2. Perform ICA
# Select most variant CpG sites
NmostVar = 200000 #recommended = 200000
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

# Extract the most contributing samples for each MC based on absolute value of contribution 
MC_active_sample = mc.activ.sample(MC_object, method = c("absolute", "reference")[1],bval = bval , MC_contrib_CpG = > MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1), output.directory = output.directory)

# Extract the most contributing samples for each MC based on differential methylation level with reference sample (here normal samples)
MC_active_sample = mc.activ.sample(MC_object, method = c("absolute", "reference")[2],bval = bval , MC_contrib_CpG = > MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1), ref = grep("N", colnames(bval), value = TRUE), output.directory = output.directory)



# Represent methylation changes in the most contributing samples vs. reference samples
mc.change(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref = grep("N", colnames(bval), value = TRUE), output.directory = output.directory)




### 4. Association of MCs with (epi)genomic characteristics
# Compute and represent enrichment of most contributing CpG sites within specific CGI-based features, chromatin states & methylation contexts.
enrich.CpG.feature(MC_object, MC_contrib_CpG, output.directory = output.directory, CpG_feature = CpG_feature)

# Compute and represent enrichment of 48 CpG category in Zhou, W., Dinh, H.Q., Ramjan, Z., Weisenberger, D.J., Nicolet, C.M., Shen, H., Laird, P.W., and Berman, B.P. (2018). DNA methylation loss in late-replicating domains is linked to mitotic cell division. Nature Genetics 50, 591–602.

# create table with catagorie
CpG_feature$nb_flanking_CpG_reccod = CpG_feature$nb_flanking_CpG
CpG_feature$nb_flanking_CpG_reccod[which(CpG_feature $nb_flanking_CpG_reccod>3)] = 3

enrich_context = matrix(NA, nrow = nrow(CpG_feature), ncol = 4*4*3)
rownames(enrich_context) = rownames(CpG_feature)
colnames(enrich_context) = paste(rep(c("HMD", "PMD", "LMR", "UMR"), each = 12), rep(rep(3:0, each = 3), length = 48), rep(c("SCGS", "SCGW", "WCGW"), length = 48), sep = "_")

i=1
for(domain_type in c("HMD", "PMD", "LMR", "UMR")){
	for(nb_flanking_CpG in c(3, 2, 1, 0)){
		for(context in c("SCGS", "SCGW", "WCGW")){
			enrich_context[which(CpG_feature$CpG_context == domain_type & CpG_feature$nb_flanking_CpG_reccod == nb_flanking_CpG & CpG_feature$context == context),i] = "yes"
			i=i+1
		}
	}
}
enrich_context = rbind(enrich_context, rep("yes", ncol(enrich_context)))
table_enrich_context = apply(enrich_context, 2, table)[colnames(enrich_context)]-1

# compute enrichment in each component
nbComp = length(MC_active_sample)

matrice_enrich = matrix(NA, ncol = nbComp, nrow = 48)
rownames(matrice_enrich) = names(table_enrich_context)
colnames(matrice_enrich) = paste0("MC", 1: nbComp)

for(i in 1:nbComp){
	context = enrich_context[MC_contrib_CpG[[i]],]
	context = rbind(context, rep("yes", ncol(context)))
		
	table_context = apply(context, 2, table)[colnames(enrich_context)]-1
	matrice_enrich[,i] = (table_context / sum(table_context)) / (table_enrich_context / sum(table_enrich_context))

}

# scale for representation
matrice_enrich = matrice_enrich-1 
matrice_enrich[which(matrice_enrich>5)] = 5
matrice_enrich[which(matrice_enrich <0)] = 0

cst_color <- colorRampPalette(c("white","white", "grey40"))(100 + 1)

pdf(file.path(output.directory, "CpG_context_Zhou.pdf"),width=8,height=8, bg = "transparent")
par(oma = c(0,0,0,0), xpd=TRUE, col = "white", mar=c(0,0,0,0), bg = "transparent")	
corrplot(matrice_enrich, col = cst_color, is.corr = FALSE, tl.col = "black", tl.cex = 1.1, mar=c(0,0,0,0), bg = "transparent",cl.lim = c(0,6))
dev.off()



### 5. Association of MCs with sample annotations
Samples_association = mc.annot(MC_object, annot = annot , save = TRUE, output.directory = output.directory)

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





