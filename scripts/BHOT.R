## THIS IS THE MASTER SCRIPT

## Load dependencies
library_list <- c("archetypes", "cowplot", "DESeq2", "dplyr", "ggplot2", "ggpubr", "ggrepel", "knitr", "MASS", "ModelMetrics", "MLeval", "nnet", "ordinal", "plyr", "predtools", "pROC", "RUVSeq", "RCRnorm", "NormqPCR", "smotefamily", "stringr")
missing_libraries <- library_list[!(library_list %in% installed.packages()[,"Package"])]

#BiocManager::install("RUVSeq")

if(length(missing_libraries) > 0) {
	cat("The following libraries are missing:\n", missing_libraries)
}

## Load all packages
lapply(library_list, library, quietly=TRUE, character.only=TRUE)

## Output data frame with all required packages
package_dump <- function() {
    si <- data.frame(sessioninfo::package_info(pkgs = c("attached")[1], include_base = FALSE, dependencies = NA))
    dump("si","histomx_pckgs.Rdmpd")
}

##------------------------------
## geometric mean
geoMean <- function(x) {
    x[x<1] <- 1;
    exp(mean(log(na.omit(x))));
}

##-------------------------------------------------
## output model performance metrics
modelEval <- function(true_values, predicted_values, threshold=NULL, print.results=TRUE) {
	
	predicted_values = predicted_values[ ! is.na(true_values) ]
	true_values = true_values[ ! is.na(true_values) ]
	true_values = true_values[ ! is.na(predicted_values) ]
	predicted_values = predicted_values[ ! is.na(predicted_values) ]
	
	roc <- pROC::roc(true_values, predicted_values, quiet=TRUE)
	
	result = list()
	result$brier_score <-round(DescTools::BrierScore(true_values, predicted_values), 3)
	result$log_loss <- round(MLmetrics::LogLoss(true_values, predicted_values) / 100, 3)
		
	pred_tab <- data.frame(truth=true_values, predicted=predicted_values)
	result$PRAUC <- round(PRROC::pr.curve(pred_tab[pred_tab$truth==1,"predicted"], pred_tab[pred_tab$truth==0,"predicted"], curve=TRUE)$auc.integral, 3)

	# the baseline in a PR curve is a horizontal line with height equal to the proportion of positive samples
	# baseline to beat is the always-positive classifier rather than any random classifier
	# Precision of a No-Skill model is equal to the fraction of positive class in the dataset
	result$PR_base <- sum(true_values) / length(true_values)
	# compare PRAUC to chance precision baseline (ie. rather than 0.5)
	#result$PRAUC / result$PRAUC_base 
	
	result$ROCAUC <- round(as.numeric(ci(roc, of="auc"))[2], 3)
	result$ROCAUC_CI_low <- round(as.numeric(ci(roc, of="auc"))[1], 3)
	result$ROCAUC_CI_high <- round(as.numeric(ci(roc, of="auc"))[3], 3)
	
	## assign class labels based on specified or Youden threshold
	if (is.null(threshold)) {
		cat("Calculating Youden threshold")
		roc_coords <- coords(roc, x="best", input="threshold", best.method="youden", transpose=F)
		predicted_outcomes <- ifelse(predicted_values > roc_coords$threshold, 1, 0)
		result$threshold <- round(roc_coords$threshold, 3)
	} else {
		predicted_outcomes <- ifelse(predicted_values > threshold, 1, 0)
		result$threshold <- threshold
	}

	# –1 = perfect misclassification and +1 = perfect classification
	result$MCC <- round(ModelMetrics::mcc(true_values, predicted_values, cutoff=result$threshold), 3)
	result$accuracy = round(sum(true_values==predicted_outcomes)/length(true_values), 3)
	
	if ( print.results ) {
		cat("Total cases that are not NA: ",length(true_values),"\n",sep="")
		cat("Accuracy=TP+TN/total: ", sum(true_values==predicted_outcomes),
		    "(",signif(result$accuracy,3),"%)\n",sep="")
	}
	
	TP = sum(true_values==1 & predicted_outcomes==1)
	TN = sum(true_values==0 & predicted_outcomes==0)
	FP = sum(true_values==0 & predicted_outcomes==1)
	FN = sum(true_values==1 & predicted_outcomes==0)
	P = TP+FN
	N = FP+TN
	precision = signif(TP/(TP+FP), 3)
	recall = signif(TP/(TP+FN), 3)
	result$PPV = signif(TP/(TP+FP), 3)
	result$TPR = signif(TP/P, 3)
	result$TNR = signif(TN/N, 3)
	result$FDR = signif(FP/(TP+FP), 3)
	result$FPR = signif(FP/N, 3)
	result$balanced_accuracy = signif((result$TPR + result$TNR)/2, 3) #balanced accuracy
	result$F1 = signif(2*((precision*recall)/(precision+recall)), 3) #harmonic mean
	
	if ( print.results ) {
		cat("PPV: (precision)=TP/(TP+FP)= ", result$PPV,"%\n",sep="")
		cat("TPR: (sensitivity)=TP/(TP+FN)= ", result$TPR,"%\n",sep="")
		cat("TNR: (specificity)=TN/(TN+FP)= ", result$TNR,"%\n",sep="")
		cat("FDR: (false discovery)=1-PPV= ", result$FDR,"%\n",sep="")
		cat("FPR: FP/N=1-TNR= ", result$FPR,"%\n",sep="")
		cat("F1: 2(PPV*TPR)/(PPV+TPR)= ", result$F1,"%\n",sep="")
	}
	#if ( print.results ) { invisible(result) }
	return(data.frame(result))
}

##-------------------------------------------------
## plot score cutoff versus TPR, TNR, Accuracy, F1
plotPrediction <- function(truth, predicted, threshold) {
	
	cutoffs <- data.frame(matrix(nrow=length(seq(0.1,1, 0.1)), ncol=5))
	colnames(cutoffs) <- c("cutoff", "TPR", "TNR", "Accuracy", "F1")
	cutoffs$cutoff <- seq(0.1,1, 0.1)
	
	for (i in 1:nrow(cutoffs)) {
		cutpoint <- cutoffs$cutoff[i]
		cutoffs[i,"TPR"] <- modelEval(truth, predicted, threshold=cutpoint, print.results=F)$TPR
		cutoffs[i,"TNR"] <- modelEval(truth, predicted, threshold=cutpoint, print.results=F)$TNR
		cutoffs[i,"Accuracy"] <- modelEval(truth, predicted, threshold=cutpoint, print.results=F)$accuracy
		cutoffs[i,"F1"] <- modelEval(truth, predicted, threshold=cutpoint, print.results=F)$F1
	}
	
	plot(cutoffs$cutoff, cutoffs$TPR, col="blue", type="l", lty=4, lwd=2, xlab="score cutoff", ylab="")
	lines(cutoffs$cutoff, cutoffs$TNR, col="red", lty=4, lwd=2)
	lines(cutoffs$cutoff, cutoffs$Accuracy, col="black", lty=3, lwd=2)
	lines(cutoffs$cutoff, cutoffs$F1, col="purple", lty=4, lwd=2)
	abline(v=threshold, col="gray", lty=4, lwd=1.5)
	legend(0.5, 0.2, legend = c("TPR", "TNR", "Accuracy", "F1"), col = c("blue", "red", "black", "purple"), 
	       lty=4, bty = "o", pt.cex=2, lwd=2, cex=0.8, text.col = "black")
}

##-------------------------------------------------
## custom caret model tuning metrics; set summaryFunction=customSummary in trainControl
customSummary <- function (data, lev = NULL, model = NULL) {
	
	pr_auc <- try(MLmetrics::PRAUC(data[, lev[2]],
				       ifelse(data$obs == lev[2], 1, 0)), silent = TRUE)
	brscore <- try(mean((data[, lev[2]] - ifelse(data$obs == lev[2], 1, 0)) ^ 2), silent = TRUE)
	rocObject <- try(pROC::roc(ifelse(data$obs == lev[2], 1, 0), data[, lev[2]],
				   direction = "<", quiet = TRUE), silent = TRUE)
	if (inherits(pr_auc, "try-error")) pr_auc <- NA
	if (inherits(brscore, "try-error")) brscore <- NA
	rocAUC <- if (inherits(rocObject, "try-error")) {
		NA
	} else {
		rocObject$auc
	}
	
	tmp <- unlist(e1071::classAgreement(table(data$obs, data$pred)))[c("diag", "kappa")]
	out <- c(Accuracy = tmp[[1]],
		 Kappa = tmp[[2]],
		 AUCROC = rocAUC,
		 AUCPR = pr_auc,
		 Brier = brscore,
		 Precision = caret:::precision.default(data = data$pred,
		 				      reference = data$obs,
		 				      relevant = lev[2]),
		 Recall = caret:::recall.default(data = data$pred,
		 				reference = data$obs,
		 				relevant = lev[2]),
		 F = caret:::F_meas.default(data = data$pred, reference = data$obs,
		 			   relevant = lev[2]))
	out
}

##------------------------------
## parse a single RCC file
importRCC <- function(rccFile) {
    
	file_name <- basename(rccFile)	
	rccFile <- readLines(rccFile)
	rccFile <- trimws(rccFile, which="both") #remove leading/trailing whitespace
	rccFile <- data.frame(line = rccFile, tag = rep("", length(rccFile)), stringsAsFactors = FALSE)
	
	## lane attributes contain Fov, BD
	tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary", "Messages")

	## assign tag to each line
	for (i in 1:length(tags)) {
        	tag_start <- which(rccFile$line == paste0("<", tags[i], ">"))
        	tag_end <- which(rccFile$line == paste0("</", tags[i], ">"))
        	rccFile$tag[tag_start:tag_end] <- tags[i]
	}

	rccFile <- rccFile[rccFile$tag != "",] ## exclude lines with empty tags
	rccFile <- rccFile[!grepl("<", rccFile$line),] ## remove lines with tag names start
	rccFile <- rccFile[!grepl("</", rccFile$line),] ## remove linear with tag names end

	## get sample ID
	samp_attr <- rccFile[rccFile$tag=="Sample_Attributes",]
	id_line <- samp_attr[grep("^ID", samp_attr$line),]
	sampID <- do.call('rbind', strsplit(as.character(id_line$line), ','))[2]

	# exclud messages; causes import issues
	all_attr <- rccFile[!rccFile$tag %in% c("Code_Summary", "Messages"),]

	parse_attr <- data.frame(do.call('rbind', strsplit(as.character(all_attr$line), ',')), stringsAsFactors=FALSE)
	parse_attr <- rbind(parse_attr, data.frame(X1="FileName", X2=file_name))
	colnames(parse_attr) <- c("variable", sampID)

	ID_lnum <- grep("^ID", parse_attr[,1])
	parse_attr[,1][ID_lnum] <- c("sampID", "laneID")
	
	## get raw counts: split columns and remove header row
	parse_counts <- rccFile[rccFile$tag=="Code_Summary",]
	parse_counts <- data.frame(do.call('rbind', strsplit(as.character(parse_counts$line), ',', fixed=TRUE)), stringsAsFactors=FALSE)
	parse_counts <- parse_counts[grep("^Code", parse_counts[,1], invert=TRUE),]
	colnames(parse_counts) <- c("CodeClass", "Name", "Accession", sampID)
	parse_counts[,sampID] <- sapply(parse_counts[,sampID], as.numeric)
    
	## return a nested list of raw counts and sample attributes
	rcc_summary <- list(counts=parse_counts, 
                        attributes=parse_attr)
	return(rcc_summary)
}

## call importRCC function to import and parse data, combine all into count table and attribute table
parseRCC <- function(rccFiles) {
    
    count_attr_list <- lapply(rccFiles, importRCC)
    
    allCounts <- do.call(rbind, count_attr_list)[,1]
    allAttributes <- do.call(rbind, count_attr_list)[,2]
    
    ## collapse raw counts to single table
    countTable <- Reduce(function(x,y) merge(x, y, all=TRUE, by=c("CodeClass", "Name", "Accession")), allCounts)
    rownames(countTable) <- countTable$Name
    attTable <- Reduce(function(x,y) merge(x, y, all=TRUE, by=c("variable")), allAttributes)
    tot_rcc <- dim(countTable)[2] - 3 ## cols 1-3=annotation columns
    
    ## Reporter Library File is used during image processing to assign target identities to barcodes, based on CodeSet name
    cat('----------------------------------', '\n');
    cat("CodeSet Name:", unique(as.character(attTable[attTable$variable=="GeneRLF",][-1])));
    cat('\n');
    cat("Successfully imported data for", tot_rcc, "sample(s)", '\n');
    cat(nrow(countTable[countTable$CodeClass=="Endogenous",]), "endogenous genes", '\n');
    cat(nrow(countTable[countTable$CodeClass=="Housekeeping",]), "housekeeping genes", '\n');
    cat(nrow(countTable[countTable$CodeClass=="Positive",]), "positive genes", '\n');
    cat(nrow(countTable[countTable$CodeClass=="Negative",]), "negative genes", '\n');
    cat('----------------------------------\n');
    
    return(list(counts=countTable,
                attributes=attTable))
}

##------------------------------
## QC RCC attributes: output binding density plot and per sample qc table
qcAttributes <- function(attributes_table) {
    
    qcTab <- attributes_table[attributes_table$variable %in% c("BindingDensity", "FovCount", "FovCounted"),]
    
    ## Imaging/Field of View measure of % of requested fields of view successfully scanned in each cartridge lane
    pct_fov <- round(as.numeric(qcTab[qcTab$variable=="FovCounted",][-1]) / as.numeric(qcTab[qcTab$variable=="FovCount",][-1]), 4)*100
    qcTab <- rbind(qcTab, c("pct_fov", pct_fov))
    
    qcTab <- setNames(data.frame(t(qcTab[,-1]), check.names=FALSE), qcTab[,1])
    qcTab$ID <- rownames(qcTab)
    rownames(qcTab) <- NULL
    
    if (any(qcTab$pct_fov < 0.75)) {
        cat(">>Imaging flag sample(s) out of range (FoV < 0.75):\n")
        print.data.frame(qcTab[qcTab$pct_fov< 0.75,c("ID", "pct_fov")])
        cat("\n")
    } else {
        cat(">>All samples pass Imaging QC\n")
    }
    
    ## Binding density: number of barcodes/μm2.
    ## recommended range is from 0.1 to 2.25 for MAX and FLEX instruments and 0.1 to 1.8 for SPRINT instruments
    ## TODO: QC based on instrument type
    if (any(qcTab$BindingDensity < 0.1 | qcTab$BindingDensity > 2.25)) {
        cat(">>Binding Density flag: sample(s) out of range (0.1 < BD < 2.25):\n")
        print.data.frame(qcTab[qcTab$BindingDensity < 0.1 | qcTab$BindingDensity > 2.25,c("ID", "BindingDensity")])
        cat("\n")
    } else {
        cat(">>All samples pass Binding Density QC\n")
    }
    
    ## per sample binding density plot
    bd_plot <- ggplot(qcTab[,c("ID", "BindingDensity")], aes(x=ID, y=BindingDensity)) +
        geom_bar(stat='identity', show.legend = FALSE) +
        xlab("") + ylab("binding density") + 
        geom_hline(yintercept=2.25, col="dodgerblue", linetype="dashed") +
        geom_hline(yintercept=0.05, col="tomato", linetype="dashed") +
        theme(panel.grid.major=element_blank(),
              panel.background=element_blank(), 
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="black", fill=NA, linewidth=1),
              plot.title = element_text(size=12, face = "bold"),
              axis.text=element_text(size=6, angle=270, family="sans", colour="black"), 
              axis.title.x=element_text(size=8, family="sans", colour="black"), 
              axis.title.y=element_text(size=8, family="sans", colour="black"))
    
    return(list(density_plot=bd_plot, 
                qcTable=qcTab))
}

##------------------------------
## Control probe quality control checks; multiple samples
## input: rows=genes, cols=CodeClass, Name, Accession, samples...
qcCounts <- function(count_table) {
	
	cat('---------------------------\n');
	cat('| Per sample control gene QC |\n');
	cat('---------------------------\n');
	
	samp_ids <- colnames(count_table[,!colnames(count_table) %in% c("CodeClass", "Name", "Accession")])
	endo_tab <- count_table[count_table$CodeClass=="Endogenous",c("Name", samp_ids)]
	hk_tab <- count_table[count_table$CodeClass=="Housekeeping",c("Name", samp_ids)]
	rownames(count_table) <- count_table$Name
	
	cat('\nOverall Assay Efficiency:\n')
	cat("(geometric mean of positive controls > 3 SD of the geometric mean of endogenous genes)\n")
	
	geo_means <- data.frame(ID=samp_ids, stringsAsFactors=FALSE)
	geo_means$endo_geo_mean <- rep("NA", nrow(geo_means))
	geo_means$hk_geo_mean <- rep("NA", nrow(geo_means))
	geo_means$pos_geo_mean <- rep("NA", nrow(geo_means))
	geo_means$neg_geo_mean <- rep("NA", nrow(geo_means))
	
	for (i in 1:length(samp_ids) ) {
		sampID <- samp_ids[i]
		endo_counts <- count_table[count_table$CodeClass=="Endogenous",][sampID][,1]
		hk_counts <- count_table[count_table$CodeClass=="Housekeeping",][sampID][,1]
		pos_counts <- count_table[count_table$CodeClass=="Positive",][sampID][,1]
		neg_counts <- count_table[count_table$CodeClass=="Negative",][sampID][,1]
		
		## add values to table
		geo_means[geo_means$ID==sampID,"endo_geo_mean"] <- round(geoMean(endo_counts), 3) # geometric mean of endogenous genes
		geo_means[geo_means$ID==sampID,"hk_geo_mean"] <- round(geoMean(hk_counts), 3) # geometric mean of HK genes
		geo_means[geo_means$ID==sampID,"pos_geo_mean"] <- round(geoMean(pos_counts), 3) # geometric mean of positive controls
		geo_means[geo_means$ID==sampID,"neg_geo_mean"] <- round(geoMean(neg_counts), 3) # geometric mean of negative controls
	}
	
	if (any(geo_means$pos_mean < geo_means$sd3)) {
		cat("Assay Efficiency is low for the following samples:\n");
		cat(geo_means[which(geo_means$pos_mean < geo_means$endo_sd3),"ID"], sep="\n")
	} else {
		cat(">>Assay Efficiency is OK for all samples\n");
	}
	
	cat('----------------------------------\n',
	    'Positive Control Linearity:\n',
	    '(correlation between the observed counts and concentrations of Positive ERCC probes)\n')
	
	## Very low assay efficiency: If the counts of all the positive controls are very low (less than ~500 even for POS_A)
	#pos_tab[pos_tab[1,]<500]
	
	pos_tab <- count_table[count_table$CodeClass=="Positive",c("Name", samp_ids)]
	
	## POS_F is considered below the limit of detection: remove when calculating linearity
	## POS_F counts are significantly higher than the negative control counts in most cases
	pos_tab <- pos_tab[grep("POS_F", pos_tab$Name, invert=TRUE),]
	
	if (!all(grepl("[[:digit:]]", pos_tab$Name))) {
		stop("Positive controls must include concentrations in their names, eg. POS_A(128)")
	}
	
	## Positive normalization factors: arithmetic mean of the geometric means of positive controls
	posTab <- data.frame(t(pos_tab[,!colnames(pos_tab) %in% "Name"]))
	posTab <- apply(posTab, MARGIN=1, FUN=geoMean);
	pos_norm_factors <- data.frame(ID=names(posTab), pos_norm_factor=mean(posTab) / posTab)
	pos_norm_factors$pos_norm_factor <- round(pos_norm_factors$pos_norm_factor, 2)
	rownames(pos_norm_factors) <- NULL
	
	if ( any(pos_norm_factors$pos_norm_factor < 0.3) | any(pos_norm_factors$pos_norm_factor > 3) ) {
		cat(">>Positive normalization factors out of range (<0.3 or >3):\n");
		print.data.frame(pos_norm_factors[pos_norm_factors$pos_norm_factor<0.3 | pos_norm_factors$pos_norm_factor>3,])
	}
	
	## extract positive control concentration
	pos_tab$concentration <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", pos_tab$Name))
	
	## calculate positive control linearity for each sample
	positiveLinearityQC <- apply(pos_tab[,samp_ids], 2, function(x) {
		round(summary(lm(x ~ pos_tab$concentration))$r.squared, 3)
	})
	
	if (any(positiveLinearityQC < 0.95)) {
		cat(">>Linearity performance is low (<0.95) for the following sample(s):\n",
		    names(which(positiveLinearityQC < 0.95)))
	} else {
		cat("\n>>Positive control linearity [lm(counts ~ concentration)] is > 0.95 for all samples\n")
	}
	
	## add positive control linearity to geoMeans table
	pcl_tab <- data.frame(PCL=positiveLinearityQC)
	pcl_tab$ID <- rownames(pcl_tab)
	
	## plot positive control counts v. concentration
	tab <- pos_tab
	
	## supress data table warning message
	suppressWarnings({
		tab <- data.table::melt(tab, measure.vars=samp_ids, verbose=FALSE)
		colnames(tab) <- c("Name", "concentration", "ID", "count")
	})
	
	linearity_plot <- ggplot(tab, aes(x=concentration, y=count, group=ID, color=ID)) +
		geom_point(shape = 21, size=3) + 
		xlab("raw counts") + ylab("concentration (fM)") +
		geom_text_repel(data=tab, aes(x=concentration, y=count, label=Name), size=4, colour="darkgray") +
		geom_smooth(aes(fill=ID), method = "lm", fullrange=TRUE, se=TRUE, size=1, 
			    color="slategray", formula = y ~ x, linetype="dashed")
	
	cat('\n----------------------------------\n',
	    'Positive Control Limit of Detection:\n',
	    '(POS_E control probe counts > mean negative control + 2*SD)\n')
	
	pos_e = count_table[grep("POS_E", count_table$Name),]
	neg_raw <- count_table[count_table$CodeClass=="Negative",]
	
	ncgMean = apply(neg_raw[,samp_ids], 2, mean)
	ncgSD = apply(neg_raw[,samp_ids], 2, sd)
	lod = ncgMean + 2*ncgSD
	pos_e_counts = pos_e[,samp_ids]
	
	## POS_E counts should be > lod
	if (any(pos_e_counts < lod)) {
		cat(">>Postive control limit of detection is low for the following sample(s):\n",
		    names(pos_e_counts[which(pos_e_counts < lod)]), '\n')
	} else {
		cat(">>Positive control E counts are > LoD for all samples\n")
	}
	
	## % endogenous (n=758) and housekeeping (n=12) genes above LoD per sample
	pctLOD <- data.frame(ID=names(lod), LOD=lod)
	pctLOD$pct_endo_alod <- rep(NA, nrow(pctLOD)) #percent endo genes above LoD
	pctLOD$pct_hk_alod <- rep(NA, nrow(pctLOD)) #percent HK genes above LoD
	pctLOD$num_hk_alod <- rep(NA, nrow(pctLOD)) #tot HK genes above LoD
	
	for (i in 1:length(samp_ids)) {
		
		idx <- samp_ids[i]
		i_lod <- pctLOD[idx,"LOD"]
		
		# percent endogenous genes > LoD
		endo_counts <- endo_tab[,idx]
		pctLOD[idx,"pct_endo_alod"] <- round(length(endo_counts[which(endo_counts > i_lod)]) / length(endo_counts) * 100, 3)
		
		# percent housekeeping genes > LoD
		hk_counts <- hk_tab[,idx]
		pctLOD[idx,"num_hk_alod"] <- round(length(hk_counts[which(hk_counts > i_lod)]), 3)
		pctLOD[idx,"pct_hk_alod"] <- round(pctLOD[idx,"num_hk_alod"] / length(hk_counts) * 100, 3)
		
	}
	
	## combine stats
	out_tab <- merge(geo_means, pctLOD, by="ID")
	out_tab <- merge(out_tab, pcl_tab, by="ID")
	out_tab <- merge(out_tab, pos_norm_factors, by="ID")
	
	cat('----------------------------------\n',
	    'Limit of Detection Flag:\n',
	    '(Sample(s) with >1 housekeeping gene below the LOD and high % of endogenous genes < LOD (lod_flag))\n');
	
	## flag samples with >1 HK gene below LoD AND %endo_genes<LOD greater than the top quartile of the distribution of %endo_genes<LoD in samples with all HK genes above LoD
	ref_blod <- quantile(100 - out_tab[out_tab$num_hk_alod==12,"pct_endo_alod"], 0.75) #samps w/ all HK genes > LoD
	out_tab$pct_endo_blod <- 100 - out_tab$pct_endo_alod # %endogenous genes < LoD for non-ref samples
	out_tab$lod_flag <- ifelse(out_tab$pct_endo_blod > ref_blod & out_tab$num_hk_alod<11, 1, 0)
	
	if (any(out_tab$lod_flag==1)) {
		cat(">>Sample(s) with LOD flag:\n",
		    out_tab[out_tab$lod_flag==1,"ID"], '\n')
	} else {
		cat(">>No samples flagged for LOD QC\n")
	}
	cat('----------------------------------\n');
	
	cat('----------------------------------\n',
	    'RNA content Flag:\n',
	    '(HK gene raw counts v. the n [n=length(HKgenes)] most highly expressed endogenous genes per sample)\n');
	
	out_tab$rna_content_r2 <- NA
	for (i in 1:nrow(out_tab)) {
		i_id <- out_tab[i,"ID"]
		i_endo <- count_table[count_table$CodeClass=="Endogenous",i_id]
		i_hk <- count_table[count_table$CodeClass=="Housekeeping",i_id]
		raw_counts <- data.frame(endo=i_endo[order(i_endo, decreasing=TRUE)][1:length(i_hk)], hk=i_hk[order(i_hk, decreasing=TRUE)])
		corr = round(summary(lm(raw_counts$endo ~ raw_counts$hk))$r.squared, 3)
		out_tab[i,"rna_content_r2"] <- corr
	}
	rna_threshold = mean(out_tab$rna_content_r2) - 3*sd(out_tab$rna_content_r2)
	if (any(out_tab$rna_content_r2 < rna_threshold)) {
		cat(">>Sample(s) with RNA content flag:\n",
		    out_tab[out_tab$rna_content_r2 < rna_threshold,"ID"], '\n')
		out_tab$rna_content_flag <- ifelse(out_tab$rna_content_r2 < rna_threshold, 1, 0)
	} else {
		cat(">>No samples flagged for RNA content\n")
	}
	cat('----------------------------------\n');
	
	#out_tab$pcl_flag <- ifelse(out_tab$PCL < 0.95, 1, 0)
	
	return(list(qc_table=out_tab,
		    pos_control_counts=pos_tab,
		    linearity_plot=linearity_plot))
	
}

##------------------------------
## TODO:combine with qcCounts > sampleQC
## QC table and plot functions:
sampleStats <- function(tab) {
    
    ## endogenous gene stats
    samp_ids <- colnames(tab[,!colnames(tab) %in% c("CodeClass", "Name", "Accession")])
    tab <- tab[tab$CodeClass=="Endogenous",]
    tab <- tab[,!colnames(tab) %in% c("CodeClass", "Name", "Accession")]
    
    endoMissing <- apply(tab, MARGIN = 2, FUN = function(x) { sum(x <= 0, na.rm=TRUE) / length(x) *100; });
    endoMedian <- apply(tab, 2, FUN=median, na.rm=TRUE);
    endoMean <- apply(tab, 2, FUN=mean, na.rm=TRUE);
    endoSD <- apply(tab, 2, FUN=sd, na.rm=TRUE);
    
    stats_df <- data.frame(row.names=names(tab[,samp_ids]),
                              Median=endoMedian, 
                              Mean=endoMean, 
                              SD=endoSD,
                              Missing_pct=endoMissing)
    stats_df <- round(stats_df, 2)
    stats_df$ID <- rownames(stats_df)
    
    return(stats_df)
}

##------------------
## batch effects/outliers
#ns.norm$batch.effects #output from nanostringnorm

##------------------
## Housekeeping gene QC and selection
## gene_annotations includes "CodeClass", "Name", "Accession"
## group_labels is optional
hkQC <- function(raw_counts, gene_annotations, group_labels=NULL) {
    
    hk_genes <- gene_annotations[gene_annotations$CodeClass=="Housekeeping","Name"]
    hk_counts <- as.matrix(raw_counts[,colnames(raw_counts) %in% hk_genes])
    
    nc_genes <- gene_annotations[gene_annotations$CodeClass=="Negative","Name"]
    ncg_counts <- as.matrix(raw_counts[,colnames(raw_counts) %in% nc_genes])
    
    #endo_genes <- gene_annotations[gene_annotations$CodeClass=="Endogenous","Name"]
    #endo_counts <- as.matrix(raw_counts[,colnames(raw_counts) %in% endo_genes])
    
    ## HK genes below limit of detection: return percent of samples in which HK genes is < LOD
    ncgMean = apply(ncg_counts, 1, mean)
    ncgSD = apply(ncg_counts, 1, sd)
    lod = ncgMean + 2*ncgSD
    #lod = ncgMean
    lod_df <- data.frame(ID=names(lod), lod=lod)
    
    hk_tab <- data.frame(hk_counts)
    hk_tab$ID <- rownames(hk_tab)
    hk_tab <- merge(lod_df, hk_tab, by="ID")
    rownames(hk_tab) <- hk_tab$ID
    hk_tab$ID <- NULL
    
    hk_blod<-NULL
    for (i in 1:nrow(hk_tab)) {
    	df<-data.frame(hk_tab[i,-1] < hk_tab[i,"lod"])
    	hk_blod <- append(hk_blod, colnames(df)[df[1,]==TRUE])
    }
    
    hk_blod_df <- data.frame(table(hk_blod))
    hk_blod_df <- hk_blod_df[order(hk_blod_df$Freq, decreasing=TRUE),]
    hk_blod_df$pct <- round(hk_blod_df$Freq / nrow(hk_tab) * 100, 3)
    
    ## 00: calculate CV for each HK gene
    hk_gene_stats <- data.frame(gene=colnames(hk_counts))
    hk_gene_stats$mean = apply(hk_counts, 2, mean, na.rm=TRUE)
    hk_gene_stats$cv = apply(hk_counts, 2, sd, na.rm=TRUE) / apply(hk_counts, 2, mean, na.rm=TRUE)
    hk_gene_stats <- hk_gene_stats[order(hk_gene_stats$cv, decreasing=F),]

    ## 01: geNorm method, Vandesompele et al. 2002
    ## rank of genes from most to least stable (2 most stable cannot be ranked)
    ## calculate average pairwise variation (sd of the log transformed expression ratios) with all other genes to determine
    ## gene-stability value M: the average pairwise variation of a particular gene with all other control genes. 
    ## Genes with the lowest M values have the most stable expression
    
    hk_gn <- NormqPCR::selectHKs(hk_counts, method = "geNorm", Symbols=colnames(hk_counts), na.rm=TRUE, minNrHK=2, log=FALSE)
    
    ## rank=1 indicates most stable
    hk_rank <- data.frame(rank=names(hk_gn$rank), gene=hk_gn$rank, stringsAsFactors=FALSE)
    
    hk_m <- data.frame(n_genes=names(hk_gn$meanM), meanM=hk_gn$meanM, stringsAsFactors=FALSE) #stability values
    hk_m <- hk_m[order(hk_m$meanM, decreasing=TRUE),]
    hk_m$n_genes <- factor(hk_m$n_genes, levels=hk_m$n_genes)
    
    hk_m <- merge(hk_m, hk_rank, by.x="n_genes", by.y="rank", all=T)
    hk_m[hk_m$n_genes==2,"gene"]<-paste(hk_m[hk_m$n_genes=="1", "gene"], collapse="|")
    hk_m <- hk_m[hk_m$n_genes!="1",]
    
    ## rank least to most stable
    gp_hk <- ggplot(data=hk_m, aes(x=n_genes, y=meanM)) + geom_point(col="dodgerblue", size=4, pch=19) + 
        ggtitle("Mean HK gene stability") +  xlab("# HK genes") +
    	geom_text_repel(data=hk_m, aes(x=n_genes, y=meanM, label=gene), colour="gray", size=4) +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_line(color="gray"),
              panel.background=element_blank(), axis.line=element_line(color="white"),
              legend.title=element_blank(),
              axis.text.y=element_text(face="bold", size=12),
              axis.text.x=element_text(face="bold", size=12),
              axis.title.x=element_text(face="bold", size=12, angle=0), 
              axis.title.y=element_text(face="bold", size=12),  
              panel.border=element_rect(colour="black", fill=NA, linewidth=1))
    
    ## select HK genes to be used in normalization
    ## manually or automatically (M<0.7)
    #num_hk_genes=as.numeric(as.character(hk_m[hk_m$meanM<0.7,][1,"n_genes"]))
    #hk_genes_keep <- hk_rank[1:num_hk_genes,"gene"]
    #hk_genes_remove <- colnames(hk_counts)[!colnames(hk_counts) %in% hk_genes_keep]
    
    ## 02: use group labels for comparison
    if (length(group_labels)>0) { 
        
    	## 02: NB regression between hk_genes ~ biological conditions/outcomes of interest
    	nb_pvals = data.frame(gene=colnames(hk_counts), pval=rep(NA, ncol(hk_counts)))
        
        for (i in 1:ncol(hk_counts)){
            hk_nb <- MASS::glm.nb(as.numeric(hk_counts[,i]) ~ as.factor(group_labels))
            nb_pvals$'pval'[i] = coef(summary(hk_nb))[2,4]
            nb_pvals <- nb_pvals[order(nb_pvals$pval, decreasing=FALSE),]
            nb_pvals$padj <- p.adjust(nb_pvals$pval, method="bonferroni")
        }
        
        ## exclude genes with padj<0.05
        #hk_genes_remove <- append(hk_genes_remove, nb_pvals[nb_pvals$padj<0.05,"gene"])
        hk_list <- list(hk_lod_stats=hk_blod_df, mval_plot=gp_hk, hk_mvals=hk_m, hk_gene_stats=hk_gene_stats, nb_pvals=nb_pvals)
        
        ## 03: NormFinder method, Andersen et al. 2004; use if groups available
        ## stepwise inclusion of more control genes until the (n + 1)th gene 
        ## has no significant contribution to the newly calculated normalization factor
        #hk_counts_t <- t(hk_counts)
        #hk_nf <- NormqPCR::selectHKs(hk_counts_t, method = "NormFinder", Symbols=rownames(hk_counts_t), group=as.factor(group_labels), minNrHK = 2, log = FALSE)
        
    } else {
        hk_list <- list(hk_lod_stats=hk_blod_df, mval_plot=gp_hk, hk_mvals=hk_m, hk_gene_stats=hk_gene_stats)
    }
    
    # if (length(hk_genes_remove)<1) {
    #     cat("No housekeeping genes to exclude based on stability")
    # }
    
    return(hk_list[lengths(hk_list) != 0])
           
}

##------------------
## Normalize data

## OPTION 1: nSolver
NSnorm <- function(eset, background_correction=FALSE, code_count=FALSE, other_norm="none", take_log=TRUE) {
    
    x <- counts(eset)
    ns_fdat <- fData(eset)
    
    ## 01 CodeCount normalization (adjust each sample based on relative value to all samples)
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/code.count.normalization.R
    if (code_count == TRUE) {
    	cat(">>Performing code count (pos control) normalization\n");
    	posTab <- t(ns_counts[rownames(ns_counts) %in% ns_fdat[ns_fdat$CodeClass=="Positive","Name"],])
    	## exclude POS_F (considered below lod)
    	posTab <- posTab[,(grep("POS_F", colnames(posTab), invert=TRUE))]
    	
    	## arithmetic mean of the geometric means of positive controls
    	pos.sample <- apply(posTab, MARGIN=1, FUN=geoMean);
    	pos.norm.factors <- data.frame(ID=names(pos.sample), mean=mean(pos.sample) / pos.sample)
    	pos.norm.factors$mean <- round(pos.norm.factors$mean, 2)
    	rownames(pos.norm.factors) <- NULL
    	
    	if ( any(pos.norm.factors$mean < 0.3) | any(pos.norm.factors$mean > 3) ) {
    		cat(">>Positive normalization factors out of range (<0.3 or >3):\n");
    		print.data.frame(pos.norm.factors[pos.norm.factors$mean<0.3 | pos.norm.factors$mean>3,])
    	}
    	
    	## multiply normalization factor by raw counts
    	x <- t(apply(x, MARGIN = 1, FUN = '*', pos.norm.factors$mean));
    } else {
    	cat("\n>>Code count normalization skipped\n")
    }

    ## 02 Background correction (subtract background from each sample)
    ## recommended for experiments in which low expressing targets common
    ## mean is least conservative; max and mean.2sd are most robust to FP
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/background.normalization.R
    
    if (background_correction == TRUE) {
        x_neg <- x[rownames(x) %in% ns_fdat[ns_fdat$CodeClass=="Negative","Name"],]
        background='mean.2sd'
        
        ## calculate max or mean of negative control probes
        background.level <- if( background != 'none') {
            if (background == 'mean.2sd') {
                mean.plus.2sd <- function(x) mean(x, na.rm=TRUE) + 2 * sd(x, na.rm=TRUE);
                background.level <- apply(x_neg, MARGIN=2, FUN=mean.plus.2sd)
            } else if (background == 'max') {
                background.level <- apply(x_neg, MARGIN=2, FUN=max, na.rm=TRUE);
            }
        }
        
        background.zscores <- data.frame(background.zscore = (background.level - mean(background.level)) / sd(background.level));
        
        if ( any(background.zscores > 3)) {
            cat("\nBackground score out of range\n")
        }
        
        ## subtract background from all genes
        x <- t(apply(x, MARGIN=1, FUN='-', background.level))
        x[x<0] <- 0 ## convert negative values to zero
        
    } else {
        cat("\n>>Background correction skipped\n")
    }
    
    ## 03 SampleContent (normalize to HK genes to account for sample or RNA content ie. pipetting fluctuations)
    ## Normalize by multiplying geometric mean of housekeeping genes by each endogenous gene
    #https://github.com/chrisrgc/NanoStringNorm/blob/master/R/sample.content.normalization.R
    ## another option is to take top endo or HK genes with lowest CV (low.cv.geo.mean)
    
    hk_genes=ns_fdat[ns_fdat$CodeClass=="Housekeeping","Name"]
    
    # calculate the normalization factor: arithmetic mean of all geometric means of HK genes
    rna.content <- apply(x[rownames(x) %in% hk_genes,], MARGIN=2, FUN=geoMean);
    hk.norm.factor <- mean(rna.content) / rna.content
    
    ## flag outlier samples
    rna.content.zscore <- data.frame(rna.zscore = (rna.content - mean(rna.content)) / sd(rna.content));
    
    if (any(abs(rna.content.zscore) > 3)) {
        cat('>>SampleContent: The following samples have sample/rna content greater than \n\t3 standard deviations from the mean:\n');
        print(signif(subset(rna.content.zscore, abs(rna.content.zscore) > 3),3));
        cat('\n');
    }
    
    ## adjust data based on normalization factors
    x <- t(apply(x, MARGIN = 1, FUN = '*', hk.norm.factor));
    
    ## TODO:apply additional normalization
    if (other_norm == "quantile") {
    	#https://github.com/cran/NanoStringNorm/blob/master/R/other.normalization.quantile.R
    } else if (other_norm == "vsn") {
    	#https://github.com/cran/NanoStringNorm/blob/master/R/other.normalization.vsn.R
    	x.fit <- vsn::vsn2(x, verbose = verbose, ...);
    	x.predict <- predict(x.fit, newdata = x.predict, useDataInFit = TRUE);
    	x.predict <- data.frame(x.predict);
    }
    
    ## endogenous genes only
    #x.norm <- x.norm[rownames(x.norm) %in% ns_fdat[ns_fdat$CodeClass=="Endogenous","Name"]]
    
    if (take_log == TRUE) {
    	cat('>>Counts are log transformed')
    	x <- log2(x);
    }
    
    return(x)
}

## OPTION 2: RUV
## default uses housekeeping genes as negative controls
RUVnorm <- function(eset, k, method="RUVg", control_genes=NULL) {
    
    cat("k =", k, '\n');
    ## upper quantile normalization (Bullard 2010)
    eset <- betweenLaneNormalization(eset, which="upper")
    fdat <- fData(eset)
    endo_genes <- fdat[fdat$CodeClass=="Endogenous","Name"]
    hk_genes <- fdat[fdat$CodeClass=="Housekeeping","Name"]
    
    if ( isFALSE(all(colnames(fData(eset)) == c("CodeClass", "Name", "Accession"))) ) {
        stop("Feature data column names must be c('CodeClass', 'Name', 'Accession')")
    }
    
    ##------------------------------
    ## negative control genes = genes assumed not to be DE with respect to the covariate of interest
    if (!is.null(control_genes)){
    	cat("using custom genes as negative controls\n")
    	c_idx <- rownames(eset)[rownames(eset) %in% control_genes]
    } else {
    	cat("using housekeeping genes as negative controls\n")
    	control_genes <- hk_genes
    	c_idx <- rownames(eset)[rownames(eset) %in% control_genes]	
    }

    ## exclude non-endogenous and non-control genes
    counts_filt <- counts(eset)[rownames(eset) %in% c(endo_genes, control_genes),]
    
    fdat_filt <- AnnotatedDataFrame(data=fdat[fdat$Name %in% c(endo_genes, control_genes),])
    eset_filt <- newSeqExpressionSet(counts_filt, phenoData=pData(eset), featureData=fdat_filt)
    
    if(method=="RUVg") {
    	
        ## RUVg using negative control genes (ie. genes assumed not to be DE with respect to the covariate of interest)
        cat("applying RUVg normalization\n")
    	eset_filt <- RUVg(eset_filt, c_idx, k=k, isLog=FALSE, round=TRUE, center=TRUE, epsilon=1)
    	
    } else if(method=="RUVs") {
    	
        ## RUVs: negative control samples for which covariates of interest are constant (eg. centered counts for technical replicates)
        #cat("applying RUVs normalization\n")
        #eset <- RUVs(eset, c_idx, k=k, isLog=FALSE, round=TRUE)
        cat("RUVs is not yet implemented")
    	
    } else if(method=="RUVr") {
    	
        ## RUVr: uses residuals from a first pass GLM regression of the unnormalized counts on covariates of interest
        #my.glm <- glm(assay(eset) ~ as.factor(group_labels)))
        #glm_res <- residuals(my.glm)
        #eset <- RUVs(eset, c_idx, k=k, residuals=glm_res, isLog=FALSE, round=TRUE)
        cat("RUVr is not yet implemented")
    	
    } else {
        cat("Possible methods include c('RUVg', 'RUVs', 'RUVr')")
    }
    
    ##------------------------------
    ## output endogenous genes only
    endo_counts <- counts(eset_filt)[rownames(eset_filt) %in% endo_genes,]
    dds <- DESeqDataSetFromMatrix(endo_counts, colData=pData(eset_filt), design=~1)
    #rowData(dds) <- fData(eset)
    
    ## estimate size factors
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersionsGeneEst(dds)
    
    cts <- counts(dds, normalized=TRUE)
    disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
    mcols(dds)$dispGeneEst <- disp
    dds <- estimateDispersionsFit(dds, fitType="mean")
    
    ## variance stabilizing log transformation
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    mat <- assay(vsd)
    
    ## remove unwanted variation as estimated by RUVg
    covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))), drop=FALSE])
    mat <- removeBatchEffect(mat, covariates=covars)
    assay(vsd) <- mat
    
    ## return expression set with log transformed normalized counts and unwanted variation
    return(list(eset = eset_filt, 
                vsd = vsd))
}

## OPTION 3: RCR
RCRnorm <- function(counts) {
    
    ns.anno <- ns.counts[,c(1:3)];
    rownames(ns.anno) <- ns.anno$Name
    
    ns.data <- ns.counts[,-c(1:3)];
    rownames(ns.data) <- ns.anno$Name
    
    ## input is a list of data frames for each code class
    ns.data.rcr <- list(pos_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Positive","Name"],],
                        neg_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Negative","Name"],],
                        hk_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Housekeeping","Name"],],
                        reg_dat=ns.data[rownames(ns.data) %in% ns.anno[ns.anno$CodeClass=="Endogenous","Name"],])
    
    #dput(as.numeric(gsub(".*\\((.*)\\).*", "\\1", ns.anno[ns.anno$CodeClass=="Positive","Name"])))
    
    rcr.norm <- RCRnorm::RCRnorm(ns.data.rcr, pos_conc=log10(c(128, 32, 8, 2, 0.5, 0.125)),
                        fast_method=FALSE, iter=8000, warmup=5000)
    
    return(rcr.norm)
    
}

##------------------
## Relative log expression plot (log2(expression/median across all assays))
plotRLE <- function(x_endo, is_logged=TRUE, main=NULL, xlab=TRUE) {
    
	if (!class(x_endo)[1]=="matrix") {
		x_endo <- as.matrix(x_endo)
	}
	
	#x_endo <- counts[rownames(counts) %in% ns.counts[ns.counts$CodeClass=="Endogenous","Name"],]
    
	endo_median <- apply(x_endo, 1, median) #per gene median
    
	if (is_logged == TRUE) {
        	endo_rle = x_endo - endo_median
	} else {
        	endo_rle <- log2(x_endo/endo_median)
	}
    
	suppressWarnings({
        	endo_melt <- data.table::melt(endo_rle, measure.vars=c(1:dim(endo_rle)[2]))
	})
	endo_melt <- na.omit(endo_melt)
    
	if(xlab=="TRUE") {
    		gp <- ggplot(endo_melt, aes(x=reorder(Var2, value, median, order=TRUE), y=value, fill=Var2)) + 
    			ggtitle(main) + ylab("log2(counts/median)") +
    			geom_hline(yintercept=0, linetype="dashed", color="firebrick") +
    			geom_boxplot(outlier.size=1, outlier.shape=21, outlier.fill="whitesmoke") +
    			theme(legend.position="none",
    			axis.text.x=element_text(size=6, angle=90),
    			axis.title.x=element_blank())
	} else {
    		gp <- ggplot(endo_melt, aes(x=reorder(Var2, value, median, order=TRUE), y=value, fill=Var2)) + 
    			ggtitle(main) + ylab("log2(counts/median)") +
    			geom_hline(yintercept=0, linetype="dashed", color="firebrick") +
    			geom_boxplot(outlier.size=1, outlier.shape=21, outlier.fill="whitesmoke") +
    			theme(legend.position="none",
    			axis.text.x=element_blank(),
    			axis.title.x=element_blank())	
	}

	return(gp)
}

##--------------------------------------------------------------------------------
## Differential expression analysis:
## TODO: mean-diff plot (with and w/o control genes), pval histogram

## 1 Mixture model: stats4::MLE()
## 2 simplified model: MASS::glm.nb()
## 3 loglinear model: lm()

## feature selection: 10-50 genes

## Volcano plot
plotVolcano <- function(df.fit, p_cutoff=0.05, num_genes=20, plot_title=NULL) {
    
    df.fit$gene <- rownames(df.fit)
    
    ## DESeq2 outputs NA values for adjusted p values based on independent filtering of genes that have low counts
    ## convert these to 1
    df.fit$pvalue <- ifelse(is.na(df.fit$pvalue), 1, df.fit$pvalue)
    df.fit$padj <- ifelse(is.na(df.fit$padj), 1, df.fit$padj)

    df.fit$logP <- -log10(df.fit$pvalue)
    df.fit$logPadj <- -log10(df.fit$padj)
    
    #df.fit <- plyr::mutate(df.fit, sig=ifelse(df.fit$padj < p_cutoff & abs(df.fit$log2FoldChange)>0.5, "sig", "not"))
    df.fit <- plyr::mutate(df.fit, sig=ifelse(df.fit$padj < p_cutoff, "sig", "not"))
    df.fit <- df.fit[order(df.fit$pvalue, decreasing=FALSE),]
    
    ## label top genes by pvalue or logFC
    gp_volcano <- ggplot() + ylab("-log10(pval)") + xlab("log2FC") + ggtitle(plot_title) +
        #geom_point(data=df.fit, aes( x=log2FoldChange, y=logP), colour="slategray", size=3, alpha=0.7) +
        geom_point(data=df.fit, aes(x=log2FoldChange, y=logP, color=sig), shape=19, size=4, alpha=0.7) +
        geom_text_repel(data=head(df.fit, num_genes), aes(x=log2FoldChange, y=logP, label=gene), colour="gray", size=3) +
        scale_color_manual(values=c("steelblue", "salmon")) + 
        geom_vline(xintercept=0, linetype="dashed", linewidth=0.4) +
        geom_vline(xintercept=1, linetype="dashed", linewidth=0.2) +
    	geom_vline(xintercept=-1, linewidth=0.2, linetype="dashed") +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.background=element_blank(), axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="black", fill=NA, linewidth=1),
              plot.title = element_text(size=14, face = "bold"), legend.position="none",
              axis.text=element_text(size=14, family="sans", colour="black"), 
              axis.title=element_text(size=14, family="sans", colour="black", face="bold"))
    
    df.fit[df.fit$gene %in% c("CCL4", "CXCL11", "CXCL10", "PLA1A", "GNLY", "ROBO4", "FGFBP2",
                              "SH2D1B", "CD160", "DARC", "ROBO4", "CDH5", "CDH13", "SOST"),] #AMR
    
    df.fit[df.fit$gene %in% c("ADAMDEC1", "ANKRD22", "CTLA4", "IFNG", "ICOS", "BTLA", 
                              "CD96", "LAG3", "SIRPG", "LCP2", "DUSP2", "CD8A"),] #TCMR
    
    return(gp_volcano)
    
}

##------------------
## QC a single RCC file
## if saveFiles=TRUE, the PCL plot and attribute table will be saved to the outPath
rccQC <- function(RCCfile, outPath, saveFiles=FALSE, output_id="sample_id") {
    
    if (rlang::is_empty(outPath)) {
        cat("Output path:", outPath, "\n")
        outPath=dirname(RCCfile)
    } else {
        cat("Output path:", outPath, "\n")
        outPath=outPath
    }
    
    ## import + parse RCC file
    ns.data <- parseRCC(RCCfile)
    
    if (output_id=="sample_id") {
    	newID <- colnames(ns.data$attributes)[2]
    } else if (output_id=="file_name") {
    	newID <- gsub("\\.RCC", "", ns.data$attributes[ns.data$attributes$variable=="FileName",2])
    } else {
    	newID <- output_id
    }
    
    ## remove spaces and special characters from sample ID
    newID <- gsub(" ", "_", newID)
    newID <- gsub("[<>(),;!@#$%?&/+]", "", newID)
    
    ## extract raw counts
    countTable <- ns.data$counts
    rownames(countTable) <- countTable$Name
    colnames(countTable)[4]<-newID
    
    ## extract attribute table
    attTable <- ns.data$attributes
    rownames(attTable) <- attTable$variable
    attTable$variable <- NULL
    colnames(attTable)<-newID
    attTable <- data.frame(t(attTable), check.names=FALSE)
    
    newOut <- paste0(outPath, newID, "/")
    dir.create(newOut, recursive=TRUE, showWarnings=FALSE)
    
    numeric_vars <- c("FileVersion", "laneID", "BindingDensity", "FovCount", "FovCounted")
    attTable[,numeric_vars] <- sapply(attTable[,numeric_vars], as.numeric)
    
    ## recommended threshold: >75%
    attTable$'% registered FOVs' <- round(attTable$FovCounted/attTable$FovCount*100, 2)
    
    ## geometric means
    attTable$'geoMean POS genes' <- exp(mean(log(countTable[countTable$CodeClass=="Positive",newID]+1)))
    attTable$'geoMean NEG genes' <- exp(mean(log(countTable[countTable$CodeClass=="Negative",newID]+1)))
    attTable$'geoMean ENDO genes' <- exp(mean(log(countTable[countTable$CodeClass=="Endogenous",newID]+1)))
    attTable$'geoMean HK genes' <- exp(mean(log(countTable[countTable$CodeClass=="Housekeeping",newID]+1)))
    
    ## neg/pos control gene QC
    ## POS_E counts should be > 2 sd * mean(NEG)
    pos_e = countTable[grep("POS_E", countTable$Name),]
    ncg <- countTable[countTable$CodeClass=="Negative",]
    
    ncgMean = mean(ncg[,newID])
    ncgSD = sd(ncg[,newID])
    lod = ncgMean + 2*ncgSD
    llod = ncgMean - 2*ncgSD
    pos_e_counts = pos_e[,-c(1:3)]
    
    ## Check if any HK genes below geoMean(ncg)
    hk_exp <- countTable[countTable$CodeClass=="Housekeeping",4]
    if (any(hk_exp < lod)) {
    	attTable$'HK below LoD' <- length(which(hk_exp < lod))
    	cat(">>Housekeeping gene(s) with expression below limit of detection detected.")
    } else {
    	attTable$'HK below LoD' <- 0
    }
    
    ## POS_E counts should be > lod
    ## % genes above LoD
    endoCounts <- countTable[countTable$CodeClass=="Endogenous",c("Name", newID)]
    colnames(endoCounts) <- c("gene", "counts")
    
    attTable$'POS_E counts' <- pos_e_counts
    attTable$ncgMean <- ncgMean
    attTable$ncgSD <- ncgSD
    attTable$LoD <- lod
    attTable$'% ENDO genes above LoD' <- round(length(endoCounts[which(endoCounts$counts > lod),"gene"]) / nrow(endoCounts) * 100, 2)
    
    ## signal versus noise: S/N = ratio of geometric mean of HK / LoD
    attTable$'SNratio' = attTable$'geoMean HK genes' / attTable$LoD
    
    ## positive control linearity
    pos <- countTable[countTable$CodeClass=="Positive",c("Name", newID)]
    pos <- pos[grep("POS_F", pos$Name, invert=TRUE),] ## POS_F should be < limit of detection
    colnames(pos) <- c("Name", "Count")
    
    if (!all(grepl("[[:digit:]]", pos$Name))) {
        stop("Positive controls must have parenthesized concentrations: ex POS_A(128)")
    }
    
    pos$Conc <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", pos$Name))
    pcl <- round(summary(lm(pos$Count ~ pos$Conc))$r.squared, 4)
    attTable$'PCL' <- pcl
    
    plot_pos_linearity <- ggplot(pos, aes(x=Conc, y=Count)) + ggtitle(paste(expression(R^2), " = ", pcl)) +
        	geom_point(shape = 21, colour="darkblue", size=3, fill = "dodgerblue") +
        	xlab("concentration") + ylab("raw counts") +
        	geom_text_repel(data=pos, aes(x=Conc, y=Count, label=Name), size=4, colour="darkgray") +
        	geom_smooth(method = "lm", fullrange=TRUE, se=TRUE, linewidth=1, 
                	color="slategray", formula = y ~ x, linetype="dashed")
    if (saveFiles==TRUE) {
    	ggsave(paste0(outPath, newID, "/pcl_plot_", newID, "_", Sys.Date(), ".pdf"), plot=plot_pos_linearity, device="pdf", width=7, height=7)
    	
    }

    outTable <- attTable[!colnames(attTable) %in% c("Owner", "Comments", "SystemAPF", "laneID", "ScannerID", "StagePosition", "CartridgeBarcode", "CartridgeID")]
    outTable <- outTable[,c("FileVersion", "SoftwareVersion", "Date", "GeneRLF", 
                            "BindingDensity", "FovCount", "FovCounted", "% registered FOVs", 
                            "geoMean POS genes", "geoMean NEG genes", "geoMean HK genes", "geoMean ENDO genes", 
                            "POS_E counts", "ncgMean", "ncgSD", "LoD", "% ENDO genes above LoD", 'HK below LoD',
                            "PCL", "SNratio")]
    vars_to_round <- c("BindingDensity", "FovCount", "FovCounted", "% registered FOVs", 
                	"geoMean POS genes", "geoMean NEG genes", "geoMean HK genes", "geoMean ENDO genes", 
                	"POS_E counts", "ncgMean", "ncgSD", "LoD", "% ENDO genes above LoD", 
                	"PCL", "SNratio")
    ## round all numeric values to 2 decimal points
    outTable[vars_to_round] <- format(round(outTable[vars_to_round], 2), nsmall = 2)
    outTable <- data.frame(t(outTable), check.names=FALSE, stringsAsFactors=FALSE)
    
    if (saveFiles==TRUE) {
    	write.table(outTable, file=paste0(outPath, newID, "/run_attribute_table_", newID, "_", Sys.Date(), ".txt"),
    		    quote=FALSE, sep='\t', row.names=TRUE)
    }
    
    ##-------------------------------------------
    ## output files and plots for markdown report
    return(list(qc_table=outTable,
                pl_plot=plot_pos_linearity))
    
}

##------------------
## Differential expression analysis for comparisons of interest

## DESEq2
runDESeq <- function(raw_counts, col_data, exp_design, contrast_vector=NULL, coef_name=NULL, pAdjustMethod="BH", verbose=FALSE) {
	
	deseq.res <- DESeqDataSetFromMatrix(countData = raw_counts, 
					    colData = col_data, 
					    design = formula(exp_design))
	# Wald Test:
	# Divide LFC it by its standard error, resulting in a z-statistic
	# The z-statistic is compared to a standard normal distribution
	# p-value = the probability that a z-statistic at least as extreme as the observed value would be selected at random
	deseq.res <- DESeq(deseq.res, test="Wald", quiet=verbose)
	
	if(!is.null(coef_name)) {
		cat("Extracting results the following coefficient:\n", coef_name, "\n")
		deseq.res.sig <- results(deseq.res, name=coef_name, pAdjustMethod=pAdjustMethod, tidy=TRUE)
	} else if(!is.null(contrast_vector)) {
		cat("Extracting results the following coefficients:\n", contrast_vector, "\n")
		deseq.res.sig <- results(deseq.res, contrast=contrast_vector, pAdjustMethod=pAdjustMethod, tidy=TRUE)
	} else {
		cat("Extracting results for all coefficients:\n", resultsNames(deseq.res)[-1], "\n")
		deseq.res.sig <- results(deseq.res, pAdjustMethod=pAdjustMethod, tidy=TRUE)
	}
	
	#TODO: add option for lfcShrink() - useful for visualizing results
	
	deseq.res.sig <- deseq.res.sig[order(deseq.res.sig$pvalue, decreasing=FALSE),]
	rownames(deseq.res.sig) <- deseq.res.sig$row
	hist(deseq.res.sig$pvalue, main="")
	
	return(deseq.res.sig)
	
}


## Wilcoxon Rank-Sum test
# norm_counts: normalized count matrix (rows=genes)
# exp_design: dataframe with the columns sample ID and condition/group
# contrast: column name in exp_design indicating groups to be compared (coded as 0 or 1)

runWilcox <- function(norm_counts, exp_design, contrast) {
    
    if (!exists('contrast')) {
        print("contrast name is missing")
    }
        
    tab.wilcox <- data.frame(gene=rownames(norm_counts), check.names=FALSE)
    
    dat <- data.frame(t(norm_counts), check.names=FALSE)
    dat$condition <- exp_design[,contrast][match(rownames(dat), exp_design$ID)]
    
    tab.wilcox$pvalue <- sapply(1:nrow(tab.wilcox), function(i) {
        gene <- tab.wilcox[i,]
        dat.gene <- dat[,colnames(dat) %in% c(gene, "condition")]
        p=wilcox.test(dat.gene[,1] ~ dat.gene$condition, paired=FALSE, alternative="two.sided")$p.value

        return(p)
    })
    
    tab.wilcox$padj <- p.adjust(tab.wilcox$pvalue, method = "BH")
    tab.wilcox <- tab.wilcox[order(tab.wilcox$pval, decreasing=FALSE),]
    
    counts_1 <- dat[dat$condition=="1",]
    counts_1 <- counts_1[,!colnames(counts_1) %in% "condition"]
    counts_0 <- dat[dat$condition=="0",]
    counts_0 <- counts_0[,!colnames(counts_0) %in% "condition"]
    tab.wilcox$log2FoldChange <- colMeans(counts_1) - colMeans(counts_0)
    tab.wilcox$FC <- 2^tab.wilcox$log2FoldChange
    
    return(tab.wilcox)
    
}
