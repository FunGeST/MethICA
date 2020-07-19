library(MethICA)
library(MethICAdata) # Use devtools::install_github("FunGeST/MethICAdata") to install this package comprising test data
library(corrplot)


### 1. Load data
# Methylation matrix, CpG table & sample annotation table.
data.directory <- file.path(.libPaths(), "MethICAdata")
load(file.path(data.directory,'/data/LICAFR_methylation.Rdata'),verbose = T)

output.directory = "~/Test_MethICA/"
if(!file.exists(output.directory)){ dir.create(output.directory) }

# To create the CpG_feature file adapted to your own tissue type, see the 'feature_table_script.R' script.

### 2. Perform ICA
# Select most variant CpG sites
NmostVar = 100000 #recommended = 100000
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
MC_active_sample = mc.active.sample(MC_object, method = c("absolute", "reference")[1],bval = bval , MC_contrib_CpG = MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1))

# Extract the most contributing samples for each MC based on differential methylation level with reference sample (here normal samples)
MC_active_sample = mc.active.sample(MC_object, method = c("absolute", "reference")[2],bval = bval , MC_contrib_CpG = MC_contrib_CpG, number = round(nrow(MC_object$Sample_contrib)*0.1), ref = grep("N", colnames(bval), value = TRUE))



# Represent methylation changes in the most contributing samples vs. reference samples
mc.change(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref = grep("N", colnames(bval), value = TRUE), output.directory = output.directory)




### 4. Association of MCs with (epi)genomic characteristics
# Compute and represent enrichment of most contributing CpG sites within specific CGI-based features, chromatin states & methylation contexts.
enrich.CpG.feature(MC_object, MC_contrib_CpG, output.directory = output.directory, CpG_feature = CpG_feature)

# Compute and represent enrichment of 48 CpG categories as in Zhou W et al. (Nat Genet 2018)
# create table with categories

CpG_feature = enrich.CpG.domain(CpG_feature = CpG_feature, MC_contrib_CpG = MC_contrib_CpG, MC_active_sample = MC_active_sample)




### 5. Association of MCs with sample annotations
sample.assoc = mc.annot(MC_object, annot = annot , save = TRUE, output.directory = output.directory)

#Examples of representations for sample associations
#boxplot
boxplot(MC_object$Sample_contrib[,"MC9"]~ annot[,"CTNNB1.alt"], col = c("grey30", "grey95"), ylab = "Sample contribution", xlab = "CTNNB1 status", main = "MC9 vs CTNNB1 status")

#corrplot for univariate analysis
association.corrplot(pvaltab_uni = as.matrix(factoall(sample.assoc$pval_uni[,2:ncol(sample.assoc$pval_uni)])), pvaltab_multi = as.matrix(factoall(sample.assoc$pval_multi)))








