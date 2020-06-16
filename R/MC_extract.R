# Rd
# description >> Extract Methylation components
# argument
# item >> bval >> methylation (beta-value) matrix
# item >> nb_comp >> number of components to extract
# item >> compute_stability >> if = TRUE, compute stability of components by performing n iterations and choose the most stable iteration
# item >> nb_iteration >> number of iterations required to assess the stability of components
# item >> alg.typ >> if alg.typ == "parallel" the components are extracted simultaneously (the default). if alg.typ == "deflation" the components are extracted one at a time
# item >> method >> if method == "R" then computations are done exclusively in R (default). The code allows the interested R user to see exactly what the algorithm does. if method == "C" then C code is used to perform most of the computations, which makes the algorithm run faster. During compilation the C code is linked to an optimized BLAS library if present, otherwise stand-alone BLAS routines are compiled.
# item >> maxit >> maximum number of iterations to perform
# item >> fun >> the functional form of the G function used in the approximation to neg-entropy (see ‘details’ of fastICA fonction).
# item >> tol >> a positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged
# item >> save >> if = TRUE, save output in the output.directory
# item >> output.directory >> path to save output
# value >> MC_object containing the following elements :
# - CpG_contrib >> contribution of CpGs to each component
# - Sample_contrib >> contribution of samples to each component
# - stability >> stability of each component 
# author >> Lea Meunier
# keyword >> Independent component analysis 
#` @import fastICA
# end

mc.extract <- function(bval, nb_comp = 20, compute_stability = TRUE, nb_iteration = 100, alg.typ = "parallel", method = "C", maxit = 1000, fun = c("logcosh","exp"), tol = 10^-6, save = FALSE, output.directory){
	
	# The fastICA algorithm uses a random initialization, which generates different results each time it is used. The set.seed() function is used to get the same result each time the algorithm is used 
	fixed_init = 154
	set.seed(fixed_init)
	it_1 <- fastICA::fastICA(bval, n.comp= nb_comp, alg.typ = "parallel", method = method, maxit = maxit, fun = fun, tol = tol)
	it_1$A = t(it_1$A)
	rownames(it_1$A) = colnames(bval)
	rownames(it_1$S) = rownames(bval)

	
	if(save == TRUE){
		if(missing(output.directory)){
			stop("missing output.directory")
		}else{
			if(!file.exists(output.directory)){
				dir.create(output.directory)
			}
			if(compute_stability == TRUE ){	
				compute_stability_it = file.path(output.directory, "iterations/")
				if(!file.exists(compute_stability_it)){
					dir.create(compute_stability_it)
				}
				write.table(it_1$S, file = file.path(compute_stability_it, "it_1_S.txt"), sep = "\t")
				write.table(it_1$A, file = file.path(compute_stability_it, "it_1_A.txt"), sep = "\t")
			}
		}
	}
	
	if(compute_stability == TRUE){
		
		# fix set.seed() for each iteration 		
		vect_fixed_init = seq(fixed_init, fixed_init+(nb_iteration-1)*12, length.out = nb_iteration)
		cor_matrix = matrix(0, ncol = nb_comp, nrow = nb_iteration)
		colnames(cor_matrix) = paste0("MC", 1:nb_comp)
		rownames(cor_matrix) = paste0("it_", 1:nb_iteration)
		
		for(i in 2:nb_iteration){
			
			# perform ICA with different set.seed() to compute stability of each component
			set.seed(vect_fixed_init[i])
			it_tmp = fastICA::fastICA(bval, n.comp= nb_comp, alg.typ = "parallel", method = method, maxit = maxit, fun = fun, tol = tol)
			it_tmp$A = t(it_tmp$A)
			rownames(it_tmp$A) = colnames(bval)
			rownames(it_tmp$S) = rownames(bval)
			
			# save results of iteration
			assign(paste0("it_",i), it_tmp)
			if(save == TRUE){
				write.table(get(paste0("it_",i, ""))$S, file = file.path(compute_stability_it, paste0("it_",i, "_S.txt")), sep = "\t")
				write.table(get(paste0("it_",i, ""))$A, file = file.path(compute_stability_it, paste0("it_",i, "_A.txt")), sep = "\t")
			}

			# compute correlation between iterations

			for(j in 1:(i-1)){			
					
				# compute stability of the new iteration with precedante
				tmp_cor_A = abs(cor(get(paste0("it_",i, ""))$A, get(paste0("it_",j, ""))$A))
				tmp_cor_S = abs(cor(get(paste0("it_",i, ""))$S, get(paste0("it_",j, ""))$S))
					
				# mean between correlation at CpG and Sample level
				tmp_cor_row = (apply(tmp_cor_A,1,max) + apply(tmp_cor_S,1,max))/2
				tmp_cor_col = (apply(tmp_cor_A,2,max) + apply(tmp_cor_S,2,max))/2	
				
				# component considerate correlated when cor >0.90 for the actual (row) and precedant (col) iteration
				tmp_cor_row[tmp_cor_row>=0.90] = 1
				tmp_cor_row[tmp_cor_row <0.90] = 0

				tmp_cor_col[tmp_cor_col>=0.90] = 1
				tmp_cor_col[tmp_cor_col<0.90] = 0
					
				# Incremente total count of iteration correlate for the actual (row) and precedant (col) iteration
				cor_matrix[i,] = cor_matrix[i,] + tmp_cor_row
				cor_matrix[j,] = cor_matrix[j,] + tmp_cor_col
			}
			print(paste("iteration", i))
		}

		sum_it = apply(cor_matrix,1,sum)
		names_final_it = names(which(sum_it == max(sum_it)))[1]
		it_final = get(names_final_it)
		stability = (cor_matrix[names_final_it,]/nb_iteration)*100
		
		cor_matrix = cbind(iteration = rownames(cor_matrix), cor_matrix)
		
		if(save == TRUE){
			write.table(cor_matrix, file = file.path(output.directory, "cor_it_matrix.txt"), sep = '\t', col.names = TRUE, row.names = TRUE)
		}
		
	}else{
		stability = rep(NA, nb_comp)
		it_final = it_1
	}
	
	CpG_contrib = it_final$S
	Sample_contrib = it_final$A
	colnames(CpG_contrib) = paste0("MC", 1:nb_comp)
	colnames(Sample_contrib) = paste0("MC", 1:nb_comp)

	if(save == TRUE){
		CpG_contrib_save = cbind(CpG = rownames(CpG_contrib), CpG_contrib)
		write.table(CpG_contrib_save, file = file.path(output.directory, "CpG_contrib.txt"))	
		Sample_contrib_save = cbind(Sample = rownames(Sample_contrib), Sample_contrib)	
		write.table(Sample_contrib, file = file.path(output.directory, "Sample_contrib.txt"))		
	}

	return(list("CpG_contrib" = CpG_contrib, "Sample_contrib" = Sample_contrib, "stability" = stability))
}




# Rd
# description >> Determine the most contributing CpGs for each component
# argument
# item >> MC_object >> methylation components object returned by mc.extract
# item >> method >> choose "threshold" to select CpGs with a contribution greater than the "threshold" parameter, or "number" to select a defined number of most contributing CpGs
# item >> threshold >> threshold on contribution values used to select CpGs with the "threshold" method. Set by default to the 95th percentile of CpG contribution values.
# item >> number >> number or most contributing CpGs to select with the "number" method
# value >> List of the most contributing CpG sites for each component
# author >> Lea Meunier
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end

mc.active.CpG <- function(MC_object, method = c("threshold", "number"), threshold = quantile(abs(MC_object$CpG_contrib), probs = 0.95), number = nrow(MC_object$CpG_contrib)*0.025){
	
	CpG_contrib = MC_object$CpG_contrib
	CpG_contrib_extract = NULL
	for(comp in colnames(CpG_contrib)){
		tmp_CpG_contrib = CpG_contrib[, comp]
		
		if(method == "threshold"){
			tmp_CpG = names(tmp_CpG_contrib)[which(abs(tmp_CpG_contrib)> threshold)]
		}else if(method == "number"){
			tmp_CpG = names(tmp_S_active)[order(abs(tmp_CpG_contrib))[1:number]]
		}

		CpG_contrib_extract[[comp]] <- tmp_CpG
	}
	return(CpG_contrib_extract)
}




# Rd
# description >> Determine the most contributing samples for each component
# argument
# item >> MC_object >> methylation components object returned by mc.extract
# item >> method >> choose "absolute" to select samples with the greatest contribution to the component (in absolute value), or "reference" to select samples with the strongest methylation difference with respect to reference samples over the most contributing CpG sites
# item >> bval >> methylation (beta-value) matrix
# item >> MC_contrib_CpG >> most contributing CpG sites for each component returned by the mc.active.CpG function
# item >> number >> number of samples selected (default: top 10 percent)
# item >> ref >> list of samples to be used as reference
# value >> List of the most contributing samples for each component
# author >> Lea Meunier
# keyword >> methods
# end

mc.active.sample <- function(MC_object, method = c("absolute", "reference"),bval = NULL , MC_contrib_CpG = NULL, number = round(nrow(MC_object$Sample_contrib)*0.1), ref = NULL){
	
	active_sample = NULL
	if(method == "absolute"){
		for(comp in colnames(MC_object$Sample_contrib)){
			tmp_S_active = MC_object$Sample_contrib[, comp]
			tmp_Sample = names(tmp_S_active)[order(abs(tmp_S_active))[1:number]]

			active_sample[[comp]] <- tmp_Sample
		}
	}else if (method == "reference"){
	
	
		for(comp in colnames(MC_object$CpG_contrib)){
			tmp_CpG_active = MC_contrib_CpG[[comp]]	
		
			bval_ref = bval[tmp_CpG_active, ref]
			bval_tmp = bval[tmp_CpG_active, setdiff(colnames(bval),ref)]
			bval_tmp = bval_tmp - apply(bval_ref,1,mean)
			bval_tmp = abs(bval_tmp)
			bval_tmp = apply(bval_tmp,2,mean)
			bval_tmp_order = bval_tmp[order(bval_tmp, decreasing = T)]	

			tmp_Sample = names(bval_tmp_order[1:number])
			active_sample[[comp]] <- tmp_Sample
	}
	}
	return(active_sample)
}




# Rd
# description >> Generate graphs showing methylation changes between reference samples and the most contributing samples to each component
# argument
# item >> MC_object >> methylation components object returned by mc.extract
# item >> MC_active_sample >> most contributing samples for each component returned by the mc.active.CpG function
# item >> MC_contrib_CpG >> most contributing CpG sites for each component returned by the mc.active.sample function
# item >> bval >> methylation (beta-value) matrix
# item >> ref >> list of samples to be used as reference
# item >> output.directory >> path to save output
# author >> Lea Meunier
# keyword >> visualisation
#` @import ggplot2
#` @import cowplot
# end

mc.change <- function(MC_object, MC_active_sample, MC_contrib_CpG, bval, ref, output.directory){
	
	if(!file.exists(output.directory)){
  		dir.create(output.directory)
	}
	# data frame use for axis representation in ggplot fonction
	df <- data.frame(x1 = 0, x2 = 0, y1 = 0, y2 = 1, x3 = 1, y3 = 1, x4 = 1, y4 =0)	
	for(comp in colnames(MC_object$CpG_contrib)){	
					
		output.directory_tmp = file.path(output.directory, comp)
		if(!file.exists(output.directory_tmp)){
			dir.create(output.directory_tmp)
		}
		
		samples_active = MC_active_sample[[comp]]
		
		# Separation of samples with a positive and negative contribution that do not carry the same change in methylation.
		samples_active_plus = names(which(MC_object$Sample_contrib[samples_active,comp]>mean(MC_object$Sample_contrib[ref,comp])))
		samples_active_moins = names(which(MC_object$Sample_contrib[samples_active,comp]<mean(MC_object$Sample_contrib[ref,comp])))
	
		# Average methylation calculation for samples with a positive and negative contribution separately
		if(length(samples_active_plus)>5){
			bval_T_pos = apply(bval[MC_contrib_CpG[[comp]], samples_active_plus], 1, mean)
			df_pos = data.frame(bval_T = bval_T_pos, bval_N = apply(bval[MC_contrib_CpG[[comp]], ref], 1, mean))
			
			p1 = ggplot2::ggplot(data= df_pos) + ggplot2::ggtitle(comp) + 
			ggplot2::xlab("Methylation in reference samples") + ggplot2::ylab("Methylation in most contributing samples") +
			ggplot2::theme_classic(base_size = 15) + ggplot2::theme(legend.position='none') + ggplot2::theme(panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank()) + 
    		ggplot2::stat_binhex(bins = 25, ggplot2::aes(x= bval_N, y= bval_T, alpha=..count..), fill="darkred") +
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x2, yend=y2))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x4, yend=y4))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x3, yend=y3))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x2, y=y2, xend=x3, yend=y3))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x3, y=y3, xend=x4, yend=y4))
			
		}else{
			bval_T_pos = NA
		}
	
		if(length(samples_active_moins)>5){
			bval_T_neg = apply(bval[MC_contrib_CpG[[comp]], samples_active_moins], 1, mean)
			df_neg = data.frame(bval_T = bval_T_neg, bval_N = apply(bval[MC_contrib_CpG[[comp]], ref], 1, mean))
			
			p2 = ggplot2::ggplot(data= df_neg) + ggplot2::ggtitle(comp) + 
			ggplot2::xlab("Methylation in reference samples") + ggplot2::ylab("Methylation in most contributing samples") +
			ggplot2::theme_classic(base_size = 15) + ggplot2::theme(legend.position='none') + ggplot2::theme(panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank()) + 
    		ggplot2::stat_binhex(bins = 25, ggplot2::aes(x= bval_N, y= bval_T, alpha=..count..), fill="darkblue") +
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x2, yend=y2))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x4, yend=y4))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x1, y=y1, xend=x3, yend=y3))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x2, y=y2, xend=x3, yend=y3))+
     		ggplot2::geom_segment(data= df, mapping=ggplot2::aes(x=x3, y=y3, xend=x4, yend=y4))
		}else{
			bval_T_neg = NA
		}
	
		# concatenation of mean methylation level and corresponding ref sample methylation level
		
		if(length(samples_active_moins)>5 & length(samples_active_plus)>5){
		  pdf(file.path(output.directory_tmp, paste0("meth_change_", comp, "_positive_samples.pdf")), width = 5, height = 5)			
		   print(p1)
		  dev.off()
		  pdf(file.path(output.directory_tmp, paste0("meth_change_", comp, "_negative_samples.pdf")), width = 5, height = 5)			
		   print(p2)
		  dev.off()
		}else if(length(samples_active_plus)>5){
			pdf(file.path(output.directory_tmp, paste0("meth_change_", comp, ".pdf")), width = 5, height = 5)			
				print(p1)
			dev.off()
		}else if(length(samples_active_moins)>5){
			pdf(file.path(output.directory_tmp, paste0("meth_change_", comp, ".pdf")), width = 5, height = 5)			
				print(p2)
			dev.off()
		}

	}
	
}




# Rd
# description >> Compute and represent enrichement of the most contributing CpG features to each component within specific (epi)genomic features
# argument
# item >> MC_object >> methylation components object returned by mc.extract
# item >> MC_contrib_CpG >> most contributing CpG sites for each component returned by the mc.active.sample function
# item >> output.directory >> path to save output
# author >> Lea Meunier
# keyword >> visualisation
#` @import plotrix
# end

enrich.CpG.feature <- function(MC_object, MC_contrib_CpG, output.directory, CpG_feature){
	
	order_feature = c("island", "shore", "shelf", "out")
	color_feature =c("#7D0E3C", "#FF584E", "#FFA576", "#B1FADA")
	color_decile = c(RColorBrewer::brewer.pal(9,"Blues"), "#05224C")
	order_active = c("active", "inactive")
	color_active = c("#1B1E26", "#E7E8E2")
	ch.states <- c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies")
	ch.states.col <- c("darkred", "red"," red1", "red2", "forestgreen", "darkgreen", "darkolivegreen1", "darkolivegreen2", "orange", "orange1", "yellow2", "turquoise", "lightslateblue", "palevioletred4", "olivedrab3" ,"grey50", "grey80", "grey95")
	names(ch.states.col) <- ch.states
	order_gene_feature = c("TSS500", "body", "out")
	color_gene_feature = c("#0C255C", "#90D1F0", "#D981AC")
	order_domain_type = c("HMD", "PMD", "LMR", "UMR")
	color_domain_type = c("darkred", "indianred1", "lightblue", "royalblue4")

	CpG_feature[is.na(CpG_feature$cgi_feature),"cgi_feature"] = "out"
	#CpG_feature[is.na(CpG_feature$cgi_type),"cgi_type"] = "out"


	for(comp in colnames(MC_object$CpG_contrib)){
		
		output.directory_tmp = file.path(output.directory, comp)
		if(!file.exists(output.directory_tmp)){
			dir.create(output.directory_tmp)
		}
		
		tmp_enrich = enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "gene_feature", CpG_feature = CpG_feature)
		pdf(file.path(output.directory_tmp, paste0("Gene_feature_", comp, ".pdf")), width = 6, height = 6)
		barplot(tmp_enrich[order_gene_feature], col = color_gene_feature)
		dev.off()
		
		tmp_enrich = enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "cgi_feature", CpG_feature = CpG_feature)
		pdf(file.path(output.directory_tmp, paste0("CGI_feature_", comp, ".pdf")), width = 6, height = 6)
		barplot(tmp_enrich[order_feature], col = color_feature)
		dev.off()
		
		tmp_enrich = enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "domain", CpG_feature = CpG_feature)
		pdf(file.path(output.directory_tmp, paste0("Chromatin_state_", comp, ".pdf")), width = 6, height = 6)
		barplot(tmp_enrich[order_active], col = color_active)
		dev.off()
		
		tmp_enrich = as.numeric(enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "state", CpG_feature = CpG_feature)[ch.states])
		names(tmp_enrich) = ch.states
		tmp_enrich[is.na(tmp_enrich)]=0
		pdf(file.path(output.directory_tmp, paste0("Chromatin_domain_", comp, ".pdf")), width = 6, height = 6)
		par(bg=NA)
		
		testpos<-seq(0,350,length.out = 18)
		echelle = pretty(range(tmp_enrich, na.rm = TRUE))
		diam_echelle = max(echelle)
		echelle_print_pos = seq(from = 0.5, to = diam_echelle, by = 1)
		
		names_cst18 = names(tmp_enrich)
		names_cst18[which(tmp_enrich <1.1)]=""

		oldpar<-plotrix::polar.plot(tmp_enrich,testpos,start=90,clockwise=TRUE,lwd=4,line.col=ch.states.col, show.grid.labels = 0, labels = names_cst18, label.prop = tmp_enrich[ch.states]/diam_echelle +0.15, radlab=FALSE, grid.col ="white", radial.lim = echelle, show.radial.grid = FALSE,boxed.radial = FALSE)
	
		a <- symbols(0,0,circles=1,fg="red",lwd=1.5, inches = FALSE, ylab ="", xlab="", lty=1, add=TRUE,xaxt = "n", yaxt = "n")
	
		for(j in echelle_print_pos[which(echelle_print_pos <diam_echelle)]){
			if(j != 1){
				symbols(0,0,circles=j,fg="grey",lwd=1, inches = FALSE, ylab ="", xlab="", lty=3, add=TRUE,xaxt = "n", yaxt = "n")
			}
			if(is.element(j, seq(0, 10, 1))){ text(j*0.75, j*0.75, labels = j)}
		}
		dev.off()

		tmp_enrich = enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "decile", CpG_feature = CpG_feature)
		pdf(file.path(output.directory_tmp, paste0("Rep_timing_", comp, ".pdf")), width = 6, height = 6)
		barplot(tmp_enrich[1:10], col = color_decile)
		dev.off()
		
		tmp_enrich = enrich.test(CpG_select = MC_contrib_CpG[[comp]], column = "CpG_context", CpG_feature = CpG_feature)
		pdf(file.path(output.directory_tmp, paste0("Domain_type_", comp, ".pdf")), width = 6, height = 6)
		barplot(tmp_enrich[order_domain_type], col = color_domain_type)	
		dev.off()
	}
}




# Rd
# description >> Compute associations between methylation components and sample annotations by uni and multivariate analysis
# argument
# item >> MC_object >> methylation components object returned by mc.extract
# item >> annot >> sample annotation table
# item >> selcol >> selection of annotations to be included in univariate analysis (by default all columns will be used)
# item >> save >> if TRUE, save results to the output.directory
# item >> output.directory >> path to save output
# item >> multi_theshold >> p-value threshold to include an annotation in the multivariate analysis
# value >> returns two matrices : p-values of univariate analysis and p-values of multivariate analysis
# author >> Lea Meunier
# keyword >> association
# end

mc.annot <- function(MC_object, annot, selcol = colnames(annot), save = FALSE, output.directory, multi_theshold = 0.001){
	annot = factoall(annot)
	annotS = annot[, selcol]
	contrib = MC_object$Sample_contrib
	
 	pval <- matrix(NA, nrow = ncol(annotS), ncol = ncol(contrib))
    rownames(pval) <- selcol
    colnames(pval) <- colnames(contrib)
    
    for(comp in colnames(contrib)){
		for(annot_test in selcol){
			tmp <- broom::tidy(lm(contrib[,comp] ~ annotS[, annot_test]))
			pval[annot_test, comp] <- min(tmp[grep("annotS\\[, annot_test\\]",tmp$term),"p.value"])
		}
	}
    if(save == TRUE){
    	pval_uni = pval
    	save(pval_uni, file = file.path(output.directory,"IC_vs_clinical_annot.Rdata"))
    	pval_uni = data.frame(annot = rownames(pval_uni), pval_uni)
		write.table(pval_uni,file.path(output.directory,"IC_vs_clinical_annot_univariate.txt"),sep="\t",quote=F,row.names=F)
    }
    
        
	for(sig in colnames(pval)){

		sign.fac <- colnames(annotS)[which(pval[,sig] < multi_theshold)]

		if(length(sign.fac)>0){
			formula <- paste0("contrib[,sig] ~ annotS$", sign.fac[1])
			if(length(sign.fac)>1){
				for(fac in sign.fac[2:length(sign.fac)])	formula <- paste(c(formula," + annotS$",as.character(fac)),collapse="")
			}
			tmp <- broom::tidy(lm(as.formula(formula)))
			tmp$Signature <- sig
			if(sig==colnames(pval)[1])	multiv <- tmp	else		multiv <- rbind(multiv,tmp)
		}
	}
	multiv <- multiv[which(multiv$p.value < 0.1),]
	write.table(multiv,file.path(output.directory,"IC_vs_clinical_annot_multivariate.txt"),sep="\t",quote=F,row.names=F)

	
	pval_multi = matrix(NA, nrow = ncol(annotS), ncol = ncol(contrib))
	rownames(pval_multi) = rownames(pval)
	colnames(pval_multi) = colnames(contrib)

	for(comp in colnames(contrib)){

		Mult_analysis = data.frame(multiv)
		Mult_analysis = Mult_analysis[which(Mult_analysis$Signature == comp),]
		Mult_analysis = Mult_analysis[which(Mult_analysis$term != "(Intercept)"),]
	
		for(j in 1:nrow(Mult_analysis)){
			for(k in 1:length(selcol)){
				tmp_var = grep(selcol[k],Mult_analysis[j,"term"])
				if(length(tmp_var)>0){
					pval_multi[selcol[k], comp] = Mult_analysis[j,"p.value"]
				}
			}		
		}
	}
	save(pval_multi,file = file.path(output.directory,"IC_vs_clinical_annot_multi.Rdata"))
	annot_multi = data.frame(annot = rownames(pval_multi), pval_multi)
	write.table(annot_multi, file.path(output.directory,"IC_vs_clinical_annot_multi.txt"),sep="\t",quote=F,row.names=F)

	return(list("pval_uni" = pval_uni, "pval_multi" = pval_multi))
}













