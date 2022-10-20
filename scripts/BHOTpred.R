
#setwd('~/Desktop/histomx-v1.9/bin/')

## predict probability of rejection on new biopsy
source('../scripts/BHOT.R')

##-----------------------------------------------------
## load models
##-----------------------------------------------------
modelPath='../model_data/kidney/models/'

##----------
## archetypes model
aa_model <- get(load(paste0(modelPath, 'aa_model.rda')))

##----------
## diagnosis base models

glm.model <- get(load(paste0(modelPath, 'dx_model_glm.rda')))
c5.model <- get(load(paste0(modelPath, 'dx_model_c5.rda')))
gbm.model <- get(load(paste0(modelPath, 'dx_model_gbm.rda')))
knn.model <- get(load(paste0(modelPath, 'dx_model_knn.rda')))
lda.model <- get(load(paste0(modelPath, 'dx_model_lda.rda')))
rf.model <- get(load(paste0(modelPath, 'dx_model_rf.rda')))
svm.model <- get(load(paste0(modelPath, 'dx_model_svm.rda')))
xgb.model <- get(load(paste0(modelPath, 'dx_model_xgb.rda')))

##----------
## banff lesions

## ordinal
g_score.model <- get(load(paste0(modelPath, 'g_ordinal_model.RData')))
ptc_score.model <- get(load(paste0(modelPath, 'ptc_ordinal_model.RData')))
i_score.model <- get(load(paste0(modelPath, 'i_ordinal_model.RData')))
t_score.model <- get(load(paste0(modelPath, 't_ordinal_model.RData')))

## binary
cg0_score.model <- get(load(paste0(modelPath, 'cg0_score_model.RData')))
v0_score.model <- get(load(paste0(modelPath, 'v0_score_model.RData')))
iifta0_score.model <- get(load(paste0(modelPath, 'iifta0_score_model.RData')))
ci1_score.model <- get(load(paste0(modelPath, 'ci1_score_model.RData')))
ct1_score.model <- get(load(paste0(modelPath, 'ct1_score_model.RData')))
cv1_score.model <- get(load(paste0(modelPath, 'cv1_score_model.RData')))

##----------
## import refSet Dx molecular scores + Dx + archetype cluster
dx_ref <- read.table('../model_data/kidney/tables/refset_molecular_scores_dx_aa.txt', header=TRUE, sep='\t')
rownames(dx_ref) <- dx_ref$ID

## import refSet lesion molecular scores + Banff scores
banff_ref <- read.table('../model_data/kidney/tables/refset_molecular_scores_banff.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

## import experimental design
exp_design <- read.table('../model_data/kidney/tables/exp_design.txt', header=TRUE, sep='\t')

## Define histology diagnosis categories
dx_amr <- c("Active AMR", "Chronic active AMR", "Chronic inactive AMR")
dx_tcmr <- c("Acute TCMR", "Chronic active TCMR")
dx_normal <- c("Normal or minimal changes", "Pristine")
dx_no_rejection <- unique(dx_ref$Dx)[!unique(dx_ref$Dx) %in% c(dx_amr, dx_tcmr)]

##-----------------------------------------------------
## import BHOT gene annotations
##-----------------------------------------------------
bhot_annot <- read.csv('../static/BHOT_annotations_v1.csv', check.names=FALSE, header=TRUE)
## exclude control genes (n=12)
rm_genes <- bhot_annot[bhot_annot$`Internal Reference Gene`=="+","Gene"]
bhot_annot <- bhot_annot[!bhot_annot$Gene %in% rm_genes,]

endats <- read.table('../static/ENDAT_genes.txt', check.names=FALSE, header=FALSE)
endats <- endats[endats$V1 %in% bhot_annot$Gene,]

bhot_cell_types <- c("B-cells", "Macrophages", "T-cells", "NK cells")

#  "MAPK", "TNF Family Signaling", added to "Type I Interferon Signaling"
bhot_pathways <- c("B-cell Receptor Signaling", "Chemokine Signaling", "Complement System", "mTOR", "NF-kappaB Signaling",
                   "Th1 Differentiation", "Th17 Differentiation", "Th2 Differentiation", "Treg Differentiation",
                   "Type I Interferon Signaling", "Type II Interferon Signaling", "IL6 Signaling",
		   "Cytotoxicity", "T-cell Receptor Signaling", "Toll-like Receptor Signaling") #important for TCMR (NK et CD8-granzyme activation, T cell activation, allorecognition/innate immunity)

bhot_annot <- bhot_annot[,colnames(bhot_annot) %in% c("Gene", "Cell Type", bhot_pathways)]

## table of total genes per category (for enrichment analysis)
bhot_pathway_totals <- data.frame(Pathway=apply(bhot_annot[,bhot_pathways], 2, function(x) {sum(x=="+")} ))
colnames(bhot_pathway_totals) <- c("Total")
bhot_pathway_totals$Pathway <- rownames(bhot_pathway_totals)

bhot_cell_type_totals <- data.frame(table(bhot_annot[bhot_annot$`Cell Type` %in% bhot_cell_types,"Cell Type"]))
colnames(bhot_cell_type_totals) <- c("CellType", "Total")
bhot_cell_type_totals <- rbind(bhot_cell_type_totals,
                               data.frame(CellType="Endothelial", Total=length(endats), check.names=F))

##-------------
## reference set RCC files (for data normalization)
refRCCpath="../refRCCs/"
refRCC <- list.files(refRCCpath, pattern=".RCC", full.names=TRUE, recursive=TRUE)

#newRCC='../test_files/amr.RCC'
#outPath='~/Downloads/'
#preds<-BHOTpred(newRCC, outPath)

##-----------------------------------------------------
## generate new predictions for a single RCC file
##-----------------------------------------------------
BHOTpred <- function(newRCC, outPath, saveFiles=FALSE) {

    cat(">>Generating HistoMx molecular report...\n")
    cat("
    O---o
     O-o
      O
     o-O
    o---O
    O---o
     O-o
      O
     o-O
    o---O\n")

    ## if no output directory provided, output files to RCC file directory
    if (exists('outPath')) {
      outPath=outPath
      cat("Output path:\n", outPath, "\n")
    } else {
      outPath=paste0(dirname(newRCC), "/")
      cat("Output path:\n", outPath, "\n")
    }

    ## get new sample ID
    new.ns.data <- parseRCC(newRCC)
    newID <- colnames(new.ns.data$attributes)[2]
    
    ## verify BHOT sequencing panel
    newRLF <- new.ns.data$attributes[new.ns.data$attributes$variable=="GeneRLF",newID]

    if (newRLF != "NS_Hs_Transplant_v1.0") {
      stop("The Gene RLF must be the same as the reference set: nCounter Human Organ Transplant Panel (NS_Hs_Transplant_v1.0)")
    } else {
      cat(">>Generating predictions for sample:", newID, "\n")
    }

    ## create sample directory in outPath
    dir.create(paste0(outPath, newID), recursive=TRUE, showWarnings=FALSE)
    newOut=paste0(outPath, newID)

    ##------------------
    ## combine refset and new RCC into single count table
    ## TODO: pre-parse reference data and save table > add new RCC to table
    #RCCfiles <- c(refRCC, newRCC)
    #ns.data <- parseRCC(RCCfiles)
    #countTable <- ns.data$counts
    #rownames(countTable) <- countTable$Name
    ##------------------

    ## load raw count table, parse newRCC, and merge counts
    ns.data <- read.table('../model_data/kidney/tables/refset_counts_raw.txt', sep='\t', header=TRUE, check.names=FALSE)
    ns.new <- parseRCC(newRCC)
    countTable <- merge(ns.data, ns.new$counts, by=c("CodeClass", "Name", "Accession"), all.x=TRUE)
    rownames(countTable) <- countTable$Name

    ## rename newID if exists in refset sample names: rownames(dx_ref)
    samp_ind <- grep(newID, colnames(countTable))
    if (length(samp_ind) > 1) {
    	colnames(countTable)[samp_ind][1]<-newID
    	colnames(countTable)[samp_ind][2]<-paste0(newID, "-new")
    	newID<-paste0(newID, "-new")
    }

    ##------------------------
    ## Normalize refset with new sample
    
    ## expression matrix
    ns.counts <- as.matrix(countTable[,-c(1:3)])
    
    ## feature data
    ns.anno <- countTable[,c(1:3)]
    rownames(ns.anno) <- ns.anno$Name
    
    ## phenotype data
    pdat <- data.frame(ID=colnames(ns.counts), check.names=FALSE)
    pdat <- merge(pdat, exp_design, by="ID", all.x=TRUE, sort=FALSE)
    rownames(pdat) <- pdat$ID
    
    ## eset
    pdat <- AnnotatedDataFrame(data=pdat,)
    fdat <- AnnotatedDataFrame(data=ns.anno)
    eset <- newSeqExpressionSet(ns.counts, phenoData=pdat, featureData=fdat)
    
    cat(">>Normalizing raw count data", "\n")
    
    ## RUV normalization
    ruv_norm <- RUVnorm(eset, k=2, method="RUVg")
    ns.norm <- assay(ruv_norm$vsd) ## normalized endogenous counts
    #ns.raw <- counts(ruv_norm$eset)[rownames(counts(ruv_norm$eset)) %in% endo_genes,]
    
    ## NanoString norm
    #x <- countTable[,!colnames(countTable) %in% c("CodeClass", "Name", "Accession")]
    ## multiply normalization factor by raw counts
    #x <- t(apply(x, MARGIN = 1, FUN = '*', pos.norm.factors$mean));
    ## SampleContent (normalize to HK genes to account for sample or RNA content ie. pipetting fluctuations)
    ## Normalize by substracting geometric mean of housekeeping genes from each endogenous gene
    #hk_genes=countTable[countTable$CodeClass=="Housekeeping","Name"]
    ## calculate the normalization factor: arithmetic mean of the geometric means of positive controls
    #rna.content <- apply(x[rownames(x) %in% hk_genes,], MARGIN=2, FUN=geoMean);
    #hk.norm.factor <- mean(rna.content) / rna.content
    #x.norm <- t(apply(x, MARGIN = 1, FUN = '*', hk.norm.factor));
    #x.norm <- log2(x.norm + 1);
    ## return endogenous probes
    #ns.norm <- x.norm[rownames(x.norm) %in% countTable[countTable$CodeClass=="Endogenous","Name"],]
    ##------------------------

    ## new sample(s) normalized with refSet
    new.ns.norm <- data.frame(counts=ns.norm[,newID], check.names=FALSE)
    colnames(new.ns.norm) <- newID

    ##--------------------------
    ## BKV expression: new biopsy v. normal and BKV refset samples
    ##--------------------------
    bk_genes <- c("BK  large T Ag", "BK  VP1")
    norm_bx_ids <- rownames(dx_ref[dx_ref$Dx %in% dx_normal,])
    bk_bx_ids <- rownames(dx_ref[dx_ref$Dx %in% c("BK virus nephropathy"),])

    bk_counts <- data.frame(ns.norm[,colnames(ns.norm) %in% c(norm_bx_ids, bk_bx_ids, newID)], check.names=FALSE)
    bk_counts <- bk_counts[rownames(bk_counts) %in% bk_genes,]
    rownames(bk_counts) <- gsub(".", "", rownames(bk_counts), fixed=TRUE)

    bkv_tab <- data.frame(normal_median=apply(bk_counts[,colnames(bk_counts) %in% norm_bx_ids], 1, function(x) { median(x, na.rm=TRUE) }),
                normal_sd=apply(bk_counts[,colnames(bk_counts) %in% norm_bx_ids], 1, function(x) { sd(x, na.rm=TRUE) }),
                bkv_median=apply(bk_counts[,!colnames(bk_counts) %in% bk_bx_ids], 1, function(x) { median(x, na.rm=TRUE) }),
                bkv_sd=apply(bk_counts[,colnames(bk_counts) %in% bk_bx_ids], 1, function(x) { sd(x, na.rm=TRUE) }),
                new_counts=bk_counts[,colnames(bk_counts) %in% newID])
    bkv_tab <- data.frame(apply(bkv_tab, 2, round, digits=3), check.names = FALSE)

    #new counts > ref_median +/- n*SD(normal counts)
    #bkv_tab$upper_sd <- bkv_tab$ref_median + 3*bkv_tab$ref_sd
    #bkv_tab$lower_sd <- bkv_tab$ref_median - 3*bkv_tab$ref_sd
    #bkv_tab <- bkv_tab[bkv_tab$new_counts>bkv_tab$upper_sd | bkv_tab$new_counts<bkv_tab$lower_sd,]

    ## boxplot
    dat <- data.frame(t(bk_counts), check.names=FALSE)
    dat[] <- sapply(dat[], as.numeric)
    #dat$BK_mean <- apply(dat[,1:2], 1, mean)
    
    dat$group <- ifelse(rownames(dat) %in% norm_bx_ids, "normal",
                        ifelse(rownames(dat) %in% bk_bx_ids, "BKV",
                               ifelse(rownames(dat) %in% newID, "new", rownames(dat))))
    dat$group <- factor(dat$group, levels=c("normal", "BKV", "new"))
    suppressWarnings({dat.m <- reshape2::melt(dat, id.vars="group")})
    bkv_boxplot <- ggplot(dat.m, aes(x = forcats::fct_rev(group), y = value, fill=group)) + geom_boxplot() +
      xlab("") + ylab("normalized expression") +
      scale_fill_manual(values=c("navy",  "lightsteelblue", "orange")) +
      theme(legend.position="none", panel.border=element_rect(colour="gray", fill=NA, size=1),
            axis.text=element_text(size=16, color="black"), axis.title=element_text(size=16, color="black"),
            panel.background=element_blank(), panel.grid.minor=element_line(colour="gray"))

    ## Wilcoxon rank sum test
    new_exp <- c(dat[dat$group %in% c("new"),"BK  VP1"], dat[dat$group %in% c("new"),"BK  large T Ag"])
    bkv_exp <- c(dat[dat$group %in% c("BKV"),"BK  VP1"], dat[dat$group %in% c("BKV"),"BK  large T Ag"])
    norm_exp <- c(dat[dat$group %in% c("normal"),"BK  VP1"], dat[dat$group %in% c("normal"),"BK  large T Ag"])
    
    suppressWarnings({
    new_greater_norm_pval <- wilcox.test(new_exp, norm_exp, paired=FALSE, alternative="greater")$p.value
    new_greater_bkv_pval <- wilcox.test(new_exp, bkv_exp, paired=FALSE, alternative="greater")$p.value
    new_less_bkv_pval <- wilcox.test(new_exp, bkv_exp, paired=FALSE, alternative="less")$p.value
    })
    
    # interpretation (for molecualr report)
    # alternate hypothesis=new > normal
    new_greater_normal_results <- ifelse(new_greater_norm_pval > 0.05, "BKV expression is not significantly higher than normal biopsies",
    				     ifelse(new_greater_norm_pval > 0.01 & new_greater_norm_pval <= 0.05, "BKV expression is moderately signficiantly higher than normal biopsies",
    				            "BKV expression is significantly higher than normal biopsies"))
    
    # alternate hypothesis=new > BKV
    new_greater_bkv_results <- ifelse(new_greater_bkv_pval > 0.05, "BKV expression is not signficiantly higher than BKV biopsies",
    				  ifelse(new_greater_bkv_pval > 0.01 & new_greater_bkv_pval <= 0.05, "BKV expression is moderately signficiantly higher than BKV biopsies",
    				         "BKV expression is signficiantly higher than BKV biopsies"))
    
    # alternate hypothesis=new < BKV
    new_less_bkv_results <- ifelse(new_less_bkv_pval > 0.05, "BKV expression is not signficiantly lower than BKV biopsies",
    			       ifelse(new_less_bkv_pval > 0.01 & new_less_bkv_pval <= 0.05, "BKV expression is moderately signficiantly lower than BKV biopsies",
    			              "BKV expression is signficiantly lower than BKV biopsies"))
    # return pvals based on BKV expression in new sample v. BKV
    if (median(new_exp) > median(bkv_exp)){
    	bkv_tab <- data.frame(new_greater_norm_pval=formatC(new_greater_norm_pval, format="e", digits=3),
    			      new_greater_bkv_pval=formatC(new_greater_bkv_pval, format="e", digits=3))
    	bkv_tab <- data.frame(pval=t(bkv_tab), check.names=F)
    	rownames(bkv_tab) <- c("new v. normal", "new v. BKV")
    	bkv_tab$Interpretation <- c(new_greater_normal_results, new_greater_bkv_results)
    } else {
    	bkv_tab <- data.frame(new_greater_norm_pval=formatC(new_greater_norm_pval, format="e", digits=3),
    			      new_less_bkv_pval=formatC(new_less_bkv_pval, format="e", digits=3))
    	bkv_tab <- data.frame(pval=t(bkv_tab), check.names=F)
    	rownames(bkv_tab) <- c("new v. normal", "new v. BKV")
    	bkv_tab$Interpretation <- c(new_greater_normal_results, new_less_bkv_results)
    }
    
    ##--------------------------
    ## signaling pathways
    ##--------------------------
    norm_bx_ids <- rownames(dx_ref[dx_ref$Dx %in% dx_normal,])

    new_counts <- ns.norm[,colnames(ns.norm) %in% newID]
    new_counts <- data.frame(gene=rownames(ns.norm),
                            new_counts=ns.norm[,colnames(ns.norm) %in% newID])
    new_counts <- new_counts[new_counts$new_counts!=0,]

    ## calculate mean|median expression in refset normal biopsies
    ref_counts <- ns.norm[,colnames(ns.norm) %in% norm_bx_ids]
    ref_counts <- data.frame(ref_median=apply(ref_counts, 1, function(x) { median(x, na.rm=TRUE) }),
                                ref_sd=apply(ref_counts, 1, function(x) { sd(x, na.rm=TRUE) }))
    ref_counts$gene <- rownames(ref_counts)

    ## combine new and normal counts
    new_norm_counts <- merge(new_counts, ref_counts, by="gene")

    new_norm_counts$FC <- new_norm_counts$new_counts / new_norm_counts$ref_median
    new_norm_counts <- new_norm_counts[order(new_norm_counts$FC, decreasing=TRUE),]

    #new counts > ref_median +/- n*SD(normal counts)
    new_norm_counts$upper <- new_norm_counts$ref_median + 2*new_norm_counts$ref_sd
    new_norm_counts$lower <- new_norm_counts$ref_median - 2*new_norm_counts$ref_sd

    ## apply FC or SD cutoff
    #new_top_genes <- new_norm_counts[new_norm_counts$FC>1.5,]
    new_top_genes <- new_norm_counts[new_norm_counts$new_counts>new_norm_counts$upper | new_norm_counts$new_counts<new_norm_counts$lower,]

    ##----------------------
    ## get pathway and cell type for each gene of interest
    ## input = new_top_fc
    ## output = pathway_table and cell_type_table
    if (nrow(new_top_genes) == 0) {

      ## if no genes of interest output table with 0 counts
      pathway_table <- data.frame('Pathway'=bhot_pathways, Count=0, check.names=FALSE)
      cell_type_table <- data.frame(CellType=bhot_cell_types, Count=0, check.names=FALSE)

    } else {

        new_top_genes$pathways <- NA
        new_top_genes$cell_types <- NA

        positive_pathways <- NULL
        positive_cell_types <- NULL

        for (i in 1:nrow(new_top_genes)) {

          gene=new_top_genes$gene[i]
          annot <- bhot_annot[bhot_annot$Gene==gene,]
          ct <- bhot_annot[bhot_annot$Gene==gene,"Cell Type"]

          p_pathways <- colnames(annot[which(annot=="+")])
          p_cell_types <- ct[ct %in% bhot_cell_types]

          ## add positive pathways/cell types to gene table
          new_top_genes[i,"pathways"] <- paste(p_pathways, collapse="|")
          new_top_genes[i,"cell_types"] <- paste(p_cell_types, collapse="|")

          ## output all pathway and cell type matches
          positive_pathways <- append(positive_pathways, p_pathways)
          positive_cell_types <- append(positive_cell_types, p_cell_types)

        }

        ##----
        ## pathways summary
        if (length(positive_pathways) > 0) {

          positive_pathways <- data.frame(table(positive_pathways))
          colnames(positive_pathways) <- c("Pathway", "Count")

          ## if all pathways represented output counts
          if (length(unique(positive_pathways$Pathway)) == length(bhot_pathways)) {
              pathway_table <- positive_pathways
          } else { #output 0
              negative_pathways <- bhot_pathways[!bhot_pathways %in% positive_pathways$Pathway]
              negative_pathways <- data.frame(Pathway=negative_pathways, Count=0)
              pathway_table <- rbind(positive_pathways, negative_pathways)
          }
          pathway_table <- pathway_table[order(pathway_table$Count, decreasing=TRUE),]

        } else {
          pathway_table <- data.frame(Pathway=bhot_pathways, Count=0)
        }

        ## cell types summary
        if (length(positive_cell_types) > 0) {

          pos_ct_df <- data.frame(table(positive_cell_types))
          colnames(pos_ct_df) <- c("CellType", "Count")

          ## Add total endothelial associated gene count
          num_endats <- data.frame(CellType='Endothelial', Count=sum(ifelse(new_top_genes$gene %in% endats, 1, 0)), check.names=FALSE)
          pos_ct_df <- rbind(pos_ct_df, num_endats)

          ## if cell types not enriched, output empty counts for these cell types
          negative_cell_types <- bhot_cell_types[!bhot_cell_types %in% pos_ct_df$CellType]
          if (length(negative_cell_types) > 0) {
            neg_ct_df <- data.frame(CellType=negative_cell_types, Count=0, check.names=FALSE)
            cell_type_table <- rbind(pos_ct_df, neg_ct_df)
          } else {
            cell_type_table <- pos_ct_df
            cell_type_table <- cell_type_table[order(cell_type_table$Count, decreasing=TRUE),]
          }

        } else {
          ## no enriched cell types
          cell_type_table <- data.frame(CellType=c(bhot_cell_types, "Endothelial"), Count=0, check.names=FALSE)
        }

    }

    ## test whether the pathway/cell type of interest contains more sig genes
    ## compared to those outside the pathway than expected by chance
    pathway_table <- merge(pathway_table, bhot_pathway_totals, by="Pathway")
    cell_type_table <- merge(cell_type_table, bhot_cell_type_totals, by="CellType")

    ## ratio of enriched genes per category
    pathway_table$GeneRatio <- round(pathway_table$Count/pathway_table$Total*100, 3)
    cell_type_table$GeneRatio <- round(cell_type_table$Count/cell_type_table$Total*100, 3)

    ##-------------
    ## Radar plot: pathways
    pathway_radar <- ggplot(pathway_table) +
      ## add custom panel grid: up to 100 or max generatio in sample
      #geom_hline(aes(yintercept = y), data.frame(y = c(0:max(pathway_table$GeneRatio)) ), color="white") +
      geom_hline(aes(yintercept = y), data.frame(y = c(0:100) ), color="white") +
      ## add bars to represent the cumulative track lengths
      geom_col(aes(x = reorder(str_wrap(Pathway, 7), GeneRatio), y=GeneRatio, fill=GeneRatio),
               position = "dodge2", show.legend = TRUE, alpha = .9) +
      # line for GeneRatio per pathway
      geom_segment(aes(x = reorder(str_wrap(Pathway, 7), GeneRatio), y = 0,
                       xend = reorder(str_wrap(Pathway, 7), GeneRatio),
                       yend = 100), linetype = "dashed", color = "navy") +
                       #yend = max(GeneRatio)), linetype = "dashed", color = "navy") +
      ## add dots to represent the count
      geom_point(aes(x=reorder(str_wrap(Pathway, 7), Count), y = Count), size = 3, color = "navy") +
      geom_text_repel(data=pathway_table[pathway_table$Count>0,], ## add count for pathways with count>0
                      aes(x=reorder(str_wrap(Pathway, 7), Count), y = Count, label=Count), size=3, colour="navy") +
      coord_polar(clip="off") +
      ## scale y axis so bars don't start in the center
      #scale_y_continuous(limits = c(-10, max(pathway_table$Total)),expand = c(0, 0), breaks = c(0, 5, 10, 20)) +
      ## set color gradient
      #scale_fill_gradientn("Number of genes", colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195")) +
      #scale_fill_gradient(low="lightblue", high="darkblue", limits=c(0,100)) +
      #scale_fill_gradientn(colours=brewer.pal(4, "Blues"), limits=c(0,100)) +
      scale_fill_gradientn(colours=c("#EFF3FF", "#BDD7E7", "#6BAED6", "darkblue"), limits=c(0,100)) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(color="navy", size=10, face="bold"),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            #plot.background=element_rect(fill="white"),
            #panel.grid.minor=element_line(colour="gray"),
            panel.grid.major=element_line(colour="gray"),
            panel.border=element_rect(colour=NA, fill=NA, size=5),
            legend.title = element_blank(),
            legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
            legend.position="right")

    ## Radar plot: cell types
    cell_type_radar <- ggplot(cell_type_table) +
      geom_hline(aes(yintercept = y), data.frame(y = c(0:100) ), color = "white") +
      geom_col(aes(x = reorder(str_wrap(CellType, 7), GeneRatio), y = GeneRatio, fill = GeneRatio),
               position = "dodge2", show.legend = TRUE, alpha = .9) +
      geom_segment(aes(x = reorder(str_wrap(CellType, 7), GeneRatio), y = 0,
                       xend = reorder(str_wrap(CellType, 7), GeneRatio),
                       yend = max(GeneRatio)), linetype = "dashed", color = "navy") +
      # Add dots to represent the count
      geom_point(aes(x=reorder(str_wrap(CellType, 7), Count), y = Count), size = 3, color = "navy") +
      geom_text_repel(data=cell_type_table[cell_type_table$Count>0,],
                      aes(x=reorder(str_wrap(CellType, 7), Count), y = Count, label=Count), size=3, colour="navy") +
      coord_polar(clip="off") +
      # Scale y axis so bars don't start in the center
      #scale_y_continuous(limits = c(-10, max(pathway_table$Total)),expand = c(0, 0), breaks = c(0, 5, 10, 20)) +
      #scale_fill_gradientn("Number of genes", colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195"))
      scale_fill_gradientn(colours=c("#EFF3FF", "#BDD7E7", "#6BAED6", "darkblue"), limits=c(0,100)) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(color="navy", size=10, face="bold"),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            #plot.background=element_rect(fill="white"),
            #panel.grid.minor=element_line(colour="gray"),
            panel.grid.major=element_line(colour="gray"),
            panel.border=element_rect(colour=NA, fill=NA, size=5),
            legend.title = element_blank(),
            legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
            legend.position="right")


    ##-------------
    ## Test for enrichment
    # my.counts<-matrix(c(pathway_table$Count[i],
    #                     pathway_table$Total[i] - pathway_table$Count[i],
    #                     nrow(new_top_genes)-pathway_table$Count[i],
    #                     nrow(bhot_annot)-pathway_table$Total[i]-nrow(new_top_genes)+pathway_table$Count[i]), 2, 2)
    # fisher.test(my.counts, alternative='greater')$p.value

    ## hypergeometric test (same as 1-sided fisher's exact)
    ## x-1 when lower.tail=F: interpretation of the p-value is P[X > x]

    ## pathways
    pathway_table$pval <- phyper(q=pathway_table$Count - 1,
                                 m=pathway_table$Total,
                                 n=nrow(bhot_annot)-pathway_table$Total,
                                 k=nrow(new_top_genes), lower.tail=FALSE)
    pathway_table <- pathway_table[order(pathway_table$pval, decreasing=F),]
    pathway_table$qval<-1.0
    
    ## adjust p value for pathways with count>0
    pathway_table[pathway_table$Count!=0,"qval"] <- round(p.adjust(pathway_table[pathway_table$Count!=0,"pval"], method="BH"), 3)

    ## cell types
    cell_type_table$pval <- phyper(q=cell_type_table$Count - 1,
                                 m=cell_type_table$Total,
                                 n=nrow(bhot_annot)-cell_type_table$Total,
                                 k=nrow(new_top_genes), lower.tail=FALSE)
    cell_type_table <- cell_type_table[order(cell_type_table$pval, decreasing=F),]
    cell_type_table$qval<-1.0
    
    ## adjust p value for cell types with count>0
    cell_type_table[cell_type_table$Count!=0,]$qval <- round(p.adjust(cell_type_table[cell_type_table$Count!=0,"pval"], method="BH"), 3)

    ##-----------------
    ## enrichment plots
    # pathway_table$GeneRatio[pathway_table$GeneRatio==0]<-0.001 ## add small value to display non-enriched pathways
    # pathway_table$Pathway <- factor(pathway_table$Pathway, levels=pathway_table$Pathway[order(pathway_table$GeneRatio)])
    #
    # pathway_plot <- ggplot(pathway_table, aes(x=Pathway, y=GeneRatio, fill=qval)) +
    #   geom_bar(stat = "identity") + xlab("") + coord_flip() +
    #   scale_colour_gradient2(limits = c(0, 1)) +
    #   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    #         panel.background=element_blank(), axis.line=element_line(colour="black"),
    #         panel.border=element_rect(colour="black", fill=NA, size=1),
    #         plot.title = element_text(size=12, face = "bold"),
    #         axis.text=element_text(size=12, family="sans", colour="black"),
    #         axis.title.x=element_text(size=12, family="sans", colour="black"),
    #         axis.title.y=element_text(size=12, family="sans", colour="black"))
    #
    # cell_type_table$GeneRatio[cell_type_table$GeneRatio==0]<-0.001 ## add small value to display non-enriched pathways
    # cell_type_table$CellType <- factor(cell_type_table$CellType, levels=cell_type_table$CellType[order(cell_type_table$GeneRatio)])
    #
    # cell_type_plot <- ggplot(cell_type_table, aes(x=CellType, y=GeneRatio, fill=qval)) +
    #   geom_bar(stat = "identity") + xlab("") + coord_flip() +
    #   scale_colour_gradient2(limits = c(0, 1)) +
    #   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    #         panel.background=element_blank(), axis.line=element_line(colour="black"),
    #         panel.border=element_rect(colour="black", fill=NA, size=1),
    #         plot.title = element_text(size=12, face = "bold"),
    #         axis.text=element_text(size=12, family="sans", colour="black"),
    #         axis.title.x=element_text(size=12, family="sans", colour="black"),
    #         axis.title.y=element_text(size=12, family="sans", colour="black"))

    ##--------------------------
    ## Diagnostic base models
    dx.genes <- colnames(coef(glm.model$finalModel))
    dx.genes <- dx.genes[grep("Intercept", dx.genes, invert=TRUE)]
    dx.genes <- gsub("`", "", dx.genes)
    
    new.ge <- new.ns.norm
    
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% dx.genes]
    
    ## base model prediction
    glm.preds <- predict(glm.model, newdata = new.ge, type='prob')
    c5.preds <- predict(c5.model, newdata = new.ge, type='prob')
    gbm.preds <- predict(gbm.model, newdata = new.ge, type='prob')
    knn.preds <- predict(knn.model, newdata = new.ge, type='prob')
    lda.preds <- predict(lda.model, newdata = new.ge, type='prob')
    rf.preds <- predict(rf.model, newdata = new.ge, type='prob')
    svm.preds <- predict(svm.model, newdata = new.ge, type='prob')
    xgb.preds <- predict(xgb.model, newdata = new.ge, type='prob')
    
    colnames(glm.preds) <- paste0(colnames(glm.preds), ".glm")
    colnames(c5.preds) <- paste0(colnames(c5.preds), ".c5")
    colnames(gbm.preds) <- paste0(colnames(gbm.preds), ".gbm")
    colnames(knn.preds) <- paste0(colnames(knn.preds), ".knn")
    colnames(lda.preds) <- paste0(colnames(lda.preds), ".lda")
    colnames(rf.preds) <- paste0(colnames(rf.preds), ".rf")
    colnames(svm.preds) <- paste0(colnames(svm.preds), ".svm")
    colnames(xgb.preds) <- paste0(colnames(svm.preds), ".xgb")
    
    all.preds <- cbind(glm.preds, c5.preds, gbm.preds, knn.preds, lda.preds, rf.preds, svm.preds, xgb.preds)
    all.preds$ID <- rownames(all.preds)
    
    ## Calculate median score across base models + medianCIs
    ## AMR
    amr.preds <- all.preds[,grep("amr", colnames(all.preds))]
    amr.preds[] <- lapply(amr.preds[], as.numeric)
    amr_ci <- apply(amr.preds, 1, function(x) { DescTools::MedianCI(x, conf.level=0.95, method="boot") })
    amr.preds$median <- amr_ci[1,] #ensemble score
    amr.preds$lwr.ci <- amr_ci[2,]
    amr.preds$upr.ci <- amr_ci[3,]
    amr.preds$lwr.ci[amr.preds$lwr.ci<0]<-0
    amr.preds$upr.ci[amr.preds$upr.ci>1]<-1

    ## TCMR
    tcmr.preds <- all.preds[,grep("tcmr", colnames(all.preds))]
    tcmr.preds[] <- lapply(tcmr.preds[], as.numeric)
    tcmr_ci <- apply(tcmr.preds, 1, function(x) { DescTools::MedianCI(x, conf.level=0.95, method="boot") })
    tcmr.preds$median <- tcmr_ci[1,] #ensemble score
    tcmr.preds$lwr.ci <- tcmr_ci[2,]
    tcmr.preds$upr.ci <- tcmr_ci[3,]
    tcmr.preds$lwr.ci[tcmr.preds$lwr.ci<0]<-0
    tcmr.preds$upr.ci[tcmr.preds$upr.ci>1]<-1
    
    ## normal: no rejection or injury diagnosis
    normal.preds <- all.preds[,grep("normal", colnames(all.preds))]
    normal.preds[] <- lapply(normal.preds[], as.numeric)
    normal_ci <- apply(normal.preds, 1, function(x) { DescTools::MedianCI(x, conf.level=0.95, method="boot") })
    normal.preds$median <- normal_ci[1,] #ensemble score
    normal.preds$lwr.ci <- normal_ci[2,]
    normal.preds$upr.ci <- normal_ci[3,]
    normal.preds$lwr.ci[normal.preds$lwr.ci<0]<-0
    normal.preds$upr.ci[normal.preds$upr.ci>1]<-1
    
    ## other: injury without rejection
    other.preds <- all.preds[,grep("other", colnames(all.preds))]
    other.preds[] <- lapply(other.preds[], as.numeric)
    other_ci <- apply(other.preds, 1, function(x) { DescTools::MedianCI(x, conf.level=0.95, method="boot") })
    other.preds$median <- other_ci[1,] #ensemble score
    other.preds$lwr.ci <- other_ci[2,]
    other.preds$upr.ci <- other_ci[3,]
    other.preds$lwr.ci[other.preds$lwr.ci<0]<-0
    other.preds$upr.ci[other.preds$upr.ci>1]<-1
    
    ## ensemble predictions + confidence intervals
    ensemble.preds <- data.frame(ID=rownames(amr.preds),
    			     amr=amr.preds$median, amr_lwr_ci=amr.preds$lwr.ci, amr_upr_ci=amr.preds$upr.ci,
    			     tcmr=tcmr.preds$median, tcmr_lwr_ci=tcmr.preds$lwr.ci, tcmr_upr_ci=tcmr.preds$upr.ci,
    			     normal=normal.preds$median, normal_lwr_ci=normal.preds$lwr.ci, normal_upr_ci=normal.preds$upr.ci,
    			     other=other.preds$median, other_lwr_ci=other.preds$lwr.ci, other_upr_ci=other.preds$upr.ci)
    rownames(ensemble.preds) <- ensemble.preds$ID
    
    new.dx.pred <- ensemble.preds 
    
    ##--------------------------
    ## Banff lesions

    ## g score
    g.genes <- names(g_score.model$coefficients)[-c(1:3)]
    g.genes <- gsub("`", '', g.genes)

    new.ge <- new.ns.norm
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% g.genes]

    new.g.pred <- data.frame(g=stats::predict(g_score.model, newdata=new.ge, type="prob"))
    colnames(new.g.pred) <- gsub(".fit.", "", colnames(new.g.pred))
    new.g.pred$ID <- newID

    ## ptc score
    ptc.genes <- names(ptc_score.model$coefficients)[-c(1:3)]
    ptc.genes <- gsub("`", '', ptc.genes)

    new.ge <- new.ns.norm
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% ptc.genes]

    new.ptc.pred <- data.frame(ptc=stats::predict(ptc_score.model, newdata=new.ge, type="prob"))
    colnames(new.ptc.pred) <- gsub(".fit.", "", colnames(new.ptc.pred))
    new.ptc.pred$ID <- newID

    ## i score
    i.genes <- names(i_score.model$coefficients)[-c(1:3)]
    i.genes <- gsub("`", '', i.genes)

    new.ge <- new.ns.norm
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% i.genes]

    new.i.pred <- data.frame(i=stats::predict(i_score.model, newdata=new.ge, type="prob"))
    colnames(new.i.pred) <- gsub(".fit.", "", colnames(new.i.pred))
    new.i.pred$ID <- newID

    ## t score
    t.genes <- names(t_score.model$coefficients)[-c(1:3)]
    t.genes <- gsub("`", '', t.genes)

    new.ge <- new.ns.norm
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% t.genes]

    new.t.pred <- data.frame(t=stats::predict(t_score.model, newdata=new.ge, type="prob"))
    colnames(new.t.pred) <- gsub(".fit.", "", colnames(new.t.pred))
    new.t.pred$ID <- newID

    ##--------------------------
    ## modify gene names for chronic lesion models
    new.ge <- new.ns.norm
    rownames(new.ge) <- gsub("-", ".", rownames(new.ge))
    rownames(new.ge) <- gsub("/", ".", rownames(new.ge))
    rownames(new.ge) <- gsub(" ", ".", rownames(new.ge))
    
    ##--------------------------
    ## predict prob ci>1 for new sample(s)
    ci1.genes <- names(ci1_score.model$coefficients)[-1] #LR

    new.ci <- t(new.ge)
    new.ci <- new.ci[,colnames(new.ci) %in% ci1.genes]

    new.pred <- stats::predict(ci1_score.model, newdata=data.frame(t(new.ci)), type="response") #LR
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ci1.pred <- new.pred
    colnames(new.ci1.pred) <- c("ci1_binary", "ID")

    ##--------------------------
    ## predict prob ct>1 for new sample(s)
    ct1.genes <- names(ct1_score.model$coefficients)[-1] #LR

    new.ct <- t(new.ge)
    new.ct <- new.ct[,colnames(new.ct) %in% ct1.genes]

    new.pred <- stats::predict(ct1_score.model, newdata=data.frame(t(new.ct)), type="response") #LR
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ct1.pred <- new.pred
    colnames(new.ct1.pred) <- c("ct1_binary", "ID")

    ##--------------------------
    ## predict prob cv>1 for new sample(s)
    cv1.genes <- names(cv1_score.model$coefficients)[-1] #LR

    new.cv <- t(new.ge)
    new.cv <- new.cv[,colnames(new.cv) %in% cv1.genes]

    new.pred <- stats::predict(cv1_score.model, newdata=data.frame(t(new.cv)), type="response") #LR

    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.cv1.pred <- new.pred
    colnames(new.cv1.pred) <- c("cv1_binary", "ID")

    ##--------------------------
    ## predict prob v>0 for new sample(s)
    v0.genes <- names(v0_score.model$coefficients)[-1] #LR

    new.v <- t(new.ge)
    new.v <- new.v[,colnames(new.v) %in% v0.genes]

    new.pred <- stats::predict(v0_score.model, newdata=data.frame(t(new.v)), type="response") #LR

    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.v0.pred <- new.pred
    colnames(new.v0.pred) <- c("v0_binary", "ID")

    ##--------------------------
    ## predict prob cg>0 for new samples
    cg0.genes <- names(cg0_score.model$coefficients)[-1] #LR

    new.cg <- t(new.ge)
    new.cg <- new.cg[,colnames(new.cg) %in% cg0.genes]

    new.pred <- stats::predict(cg0_score.model, newdata=data.frame(t(new.cg)), type="response") #LR

    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.cg0.pred <- new.pred
    colnames(new.cg0.pred) <- c("cg0_binary", "ID")

    ##--------------------------
    ## predict prob i_IFTA>0 for new samples
    iifta0.genes <- names(iifta0_score.model$coefficients)[-1] #LR

    new.iifta <- t(new.ge)
    new.iifta <- new.iifta[,colnames(new.iifta) %in% iifta0.genes]

    new.pred <- stats::predict(iifta0_score.model, newdata=data.frame(t(new.iifta)), type="response") #LR

    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.iifta0.pred <- new.pred
    colnames(new.iifta0.pred) <- c("iifta0_binary", "ID")

    ##-------------------------------------
    ## score table
    options(scipen = 999)
    join_list <- list(new.dx.pred,
    		      new.g.pred, new.ptc.pred, new.i.pred, new.t.pred, #ordinal
                      new.cg0.pred, new.iifta0.pred, new.v0.pred,
                      new.cv1.pred, new.ci1.pred, new.ct1.pred)
    tab <- Reduce(function(...) merge(..., all=TRUE), join_list)
    rownames(tab) <- tab$ID
    
    new_scores <- tab ## save scores before calculating %
    tab <- data.frame(round(tab[-1], 3) *100)
    tab$ID <- rownames(tab)
    new_scores_pct <- tab %>% dplyr::select(ID, everything())

    ## output score table
    if(saveFiles=="TRUE") {
      write.table(new_scores_pct, file=paste0(newOut, "/molecular_score_table_", newID, "_", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
    }

    options(scipen = 0)

    ##-----------------------
    ## output IQR and median scores by Dx for reference Dx
    
    ## add binary Banff lesion scores
    ref.tab <- merge(dx_ref, banff_ref[,c("ID", "iifta0_binary", "v0_binary", "cg0_binary", "ti1_binary", "cv1_binary", 
    				      "ci1_binary", "ct1_binary")], by="ID")
    rownames(ref.tab) <- ref.tab$ID
    
    dx_scores <- c("amr", "tcmr", "normal", "other", "iifta0_binary", "v0_binary", "cg0_binary", "ti1_binary", "cv1_binary", 
    	       "ci1_binary", "ct1_binary")
    #ref.tab <- ref.tab[apply(ref.tab[,-1], 1, function(x) { all(!is.na(x)) }),] #exclude samples without all molecular scores
    
    ## IQR = Q3 - Q1
    ## Q0= min and Q5=max
    ref.amr.quantiles <- apply(ref.tab[ref.tab$Dx_simple=="amr",dx_scores], 2, quantile, na.rm=TRUE)
    ref.tcmr.quantiles <- apply(ref.tab[ref.tab$Dx_simple=="tcmr",dx_scores], 2, quantile, na.rm=TRUE)
    ref.normal.quantiles <- apply(ref.tab[ref.tab$Dx_simple=="normal",dx_scores], 2, quantile, na.rm=TRUE)
    ref.other.quantiles <- apply(ref.tab[ref.tab$Dx_simple=="other",dx_scores], 2, quantile, na.rm=TRUE)
    
    ref.amr.quantiles <- data.frame(apply(ref.amr.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.amr.quantiles$Dx <- "AMR"
    ref.amr.quantiles$Q <- rownames(ref.amr.quantiles)

    ref.tcmr.quantiles <- data.frame(apply(ref.tcmr.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.tcmr.quantiles$Dx <- "TCMR"
    ref.tcmr.quantiles$Q <- rownames(ref.tcmr.quantiles)

    ref.normal.quantiles <- data.frame(apply(ref.normal.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.normal.quantiles$Dx <- "NONREJECTION"
    ref.normal.quantiles$Q <- rownames(ref.normal.quantiles)
    
    ref.other.quantiles <- data.frame(apply(ref.other.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.other.quantiles$Dx <- "INJURY"
    ref.other.quantiles$Q <- rownames(ref.other.quantiles)

    ref.score.quantiles <- rbind(ref.amr.quantiles, ref.tcmr.quantiles, ref.normal.quantiles, ref.other.quantiles)
    
    ##----------------------------------------------
    ## boxplots of reference biopsy molecular scores
    dat <- ref.tab[ref.tab$Dx_simple=="amr",]
    dat <- dat[,dx_scores]
    #amr_order <- names(ref.amr.median[order(ref.amr.median, decreasing=F)])
    #dat <- dat[,amr_order]
    colnames(dat) <- gsub("_score", "", colnames(dat))
    suppressWarnings({
      dat.m <- reshape2::melt(dat)
    })
    dat.m$value <- dat.m$value * 100

    boxplot_amr <- ggplot(data=dat.m, aes(x=reorder(variable, value, median, order=TRUE), y=value)) +
    	ggtitle("AMR reference biopsies") +
    	ylab("molecular score (%)") + xlab("") + ylim(c(0, 100)) +
    	geom_hline(yintercept=50, linetype="dashed", color="slategray") +
    	geom_boxplot(color="steelblue4", fill="whitesmoke") +
    	theme(legend.title=element_blank(),
    	      panel.grid.minor=element_line(colour="gray"),
    	      panel.grid.major=element_blank(),
    	      panel.background=element_blank(),
    	      axis.line=element_line(colour="white"),
    	      axis.text.x=element_text(face="bold", size=14, angle=0),
    	      axis.title.x=element_text(face="bold", size=14, angle=0),
    	      axis.text.y=element_text(size=14, angle=0),
    	      axis.title.y=element_text(face="bold", size=14),
    	      panel.border=element_rect(colour="gray", fill=NA, size=1))
    #ggsave(paste0(newOut, "/amr_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_amr, device="pdf", width=8, height=6)

    ##--------------------------
    ## combine refset and new molecular scores
    ## scores to input to PCA
    pca_scores <- c("amr", "tcmr", "normal")
    #pca_scores <- c("g0", "g1", "g2", "g3", "ptc0", "ptc1", "ptc2", "ptc3",
    #                "i0", "i1", "i2", "i3", "t0", "t1", "t2", "t3")

    ## keep only Bx with all molecular scores
    mscores_pca <- ref.tab[apply(ref.tab[,pca_scores], 1, function(x) { all(!is.na(x)) }),]
	
    ## add scores for new sample
    mscores_pca <- rbind(mscores_pca[,pca_scores], new_scores[,pca_scores])

    resPCA <- FactoMineR::PCA(mscores_pca[,!colnames(mscores_pca) %in% "ID"], scale.unit=FALSE, ncp=5, graph=FALSE)

    pca.df <- as.data.frame(resPCA$ind$coord)
    pca.df$ID <- rownames(pca.df)

    pca.df <- merge(pca.df, ref.tab[,c("ID", "Dx")], by="ID", all.x=TRUE)
    pca.df$Dx <- ifelse(pca.df$ID %in% new_scores$ID, "new", pca.df$Dx)

    ## simplify Dx
    pca.df$Dx <- ifelse(pca.df$Dx %in% c(dx_amr, dx_tcmr, dx_normal, "new"), pca.df$Dx, "Injury w/o rejection")
    #pca.df$Dx <- ifelse(pca.df$Dx %in% dx_normal, "No rejection/injury Dx", pca.df$Dx)
       
    ## highlight new biopsy
    pca.df <- plyr::mutate(pca.df, ref=ifelse(pca.df$Dx!="new", "ref", "new"))

    #col_vector=c("firebrick",  "blue3", "mediumpurple", "turquoise", "orangered", "dodgerblue", "salmon",
    #	     "olivedrab2","darkgreen", "black", "palegreen", "lightgray")
    col_vector=c("firebrick",  "blue3", "salmon", "dodgerblue", "orangered", "lightgray", "black", "gray36")
    	     
    pca_new_1_2 <- ggplot() +
        scale_fill_manual(values=col_vector) +
        #geom_point(data=pca.df, aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) +
        geom_point(data=pca.df[pca.df$ref=="ref",], aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=3, alpha=0.8) +
        geom_point(data=pca.df[pca.df$ref=="new",], aes(Dim.1, Dim.2), shape=23, size=4, alpha=1, fill="orange") +
        #geom_text_repel(data=pca.df[pca.df$ref=="new",], aes(Dim.1, Dim.2, label=ID), size=2, colour="orange") +
        xlab(paste("PC1 ", round(resPCA$eig[1,2]),"% of variance",sep="")) +
        ylab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
        theme(legend.position="top", legend.title=element_blank(),
              panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="white"),
              axis.text=element_text(size=12),
              axis.text.x=element_text(face="bold", size=12, angle=0),
              axis.title.x=element_text(face="bold", size=12, angle=0),
              axis.title.y=element_text(face="bold", size=12),
              panel.border=element_rect(colour="white", fill=NA, size=5))

    pca_new_2_3 <- ggplot() +
      scale_fill_manual(values=col_vector) +
      geom_point(data=pca.df[pca.df$ref=="ref",], aes(Dim.2, Dim.3, fill=Dx), shape=21, color="gray", size=3, alpha=0.7) +
      geom_point(data=pca.df[pca.df$ref=="new",], aes(Dim.2, Dim.3), shape=23, size=4, alpha=1, fill="orange") +
      xlab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
      ylab(paste("PC3 ", round(resPCA$eig[3,2]),"% of variance",sep="")) +
      theme(legend.position="none",
            panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(colour="white"),
            axis.text=element_text(size=12),
            axis.text.x=element_text(face="bold", size=12, angle=0),
            axis.title.x=element_text(face="bold", size=12, angle=0),
            axis.title.y=element_text(face="bold", size=12),
            panel.border=element_rect(colour="white", fill=NA, size=5))

    if (saveFiles=="TRUE") {
      ggsave(paste0(newOut, "/pca_dx_pc1_pc2_", newID, "_", Sys.Date(), ".pdf"), plot=pca_new_1_2, device="pdf", width=7, height=6)
    }

    ##-----------------------------
    ## KNN: k-nearest neighbors
    ##-----------------------------
    k=25
    #knn_tab <- rbind(dx_ref[,pca_scores], new_scores[,pca_scores])
    knn_tab <- mscores_pca
    knn_tab$ID <- rownames(knn_tab)
    knn_tab <- merge(knn_tab, pca.df[,c("ID", "Dx", "ref")], by="ID")

    train_ref <- knn_tab[knn_tab$ref %in% "ref",c("amr", "tcmr", "normal", "ID", "Dx")]
    test_new <- knn_tab[knn_tab$ref %in% "new",c("amr", "tcmr", "normal")]

    ## distance measure (for continuous features) or similarity measure (for categorical features)
    knn_res <- neighbr::knn(train_set=train_ref, k=k,
                          test_set=test_new,
                          comparison_measure="euclidean",
                          categorical_target="Dx",
                          return_ranked_neighbors=k, id="ID")

    ##----------------
    ## closest neighbors ranked by distance
    nn_new <- knn_res$test_set_scores #nearest neighbors

    nn_df <- data.frame(t(nn_new[2:length(nn_new)]))
    names(nn_df) <- "ID"

    ##----------------
    ## histology based diagnosis of nearest neighbors
    #length(nn_df$ID) == k
    nn_dx <- data.frame(table(knn_tab[knn_tab$ID %in% nn_df$ID,"Dx"]))
    colnames(nn_dx) <- c("Dx", "Total")
    nn_dx <- nn_dx[nn_dx$Dx!="new",]
    sum(nn_dx$Total) == k
    nn_dx$Percent <- round(nn_dx$Total/sum(nn_dx$Total) * 100, 0)

    nn_dx <- nn_dx[order(nn_dx$Total, decreasing=TRUE),]
    rownames(nn_dx) <- nn_dx$Dx
    nn_dx$Dx <- NULL

    ##----------------
    ## banff scores of nearest neighbors
    banff_df <- merge(nn_df, banff_ref[,c("ID", "g_score", "ptc_score", "cg_score", "i_score", "t_score", "v_score")], by="ID")

    g_df <- data.frame(g=table(banff_df$g_score, useNA="no"))
    ptc_df <- data.frame(ptc=table(banff_df$ptc_score, useNA="no"))
    cg_df <- data.frame(cg=table(banff_df$cg_score, useNA="no"))
    i_df <- data.frame(i=table(banff_df$i_score, useNA="no"))
    t_df <- data.frame(t=table(banff_df$t_score, useNA="no"))
    v_df <- data.frame(v=table(banff_df$v_score, useNA="no"))

    score_df <- data.frame(score=c("0", "1", "2", "3"))
    score_df <- merge(score_df, g_df, by.x="score", by.y="g.Var1", all=TRUE)
    score_df <- merge(score_df, ptc_df, by.x="score", by.y="ptc.Var1", all=TRUE)
    score_df <- merge(score_df, cg_df, by.x="score", by.y="cg.Var1", all=TRUE)
    score_df <- merge(score_df, i_df, by.x="score", by.y="i.Var1", all=TRUE)
    score_df <- merge(score_df, t_df, by.x="score", by.y="t.Var1", all=TRUE)
    score_df <- merge(score_df, v_df, by.x="score", by.y="v.Var1", all=TRUE)

    score_df[is.na(score_df)]<-0

    rownames(score_df) <- score_df$score
    nn_banff <- data.frame(round(apply(score_df[,-1], 2, function(x) {x/sum(x)*100}), 0), check.names=F) #% of Bx with given score by lesion
    colnames(nn_banff) <- gsub(".Freq", "", colnames(nn_banff))

    ##----------------
    ## mean molecular score of nearest neighbors
    ## TODO: if not all Dx in nn_dx table
    #nn_scores <- knn_tab[knn_tab$ref=="ref",c("ID", "normal", "amr", "TCMR", "ATI")]
    #nn_scores <- nn_scores[nn_scores$ID %in% nn_df$ID,]

    #nn_means <- apply(nn_scores[,-1], 2, function(x) { mean(as.numeric(x)) })
    #nn_means <- round(nn_means * 100, 2) #convert to %
    #nn_means <- data.frame(t(nn_means))
    #rownames(nn_means) <- "mean_score_knn"

    #kable(nn_means, row.names=TRUE, align='l', caption="", format="html") %>% kableExtra::kable_styling(position="left")
    #nn_dx$mean_score_knn <- t(nn_means)
    #write.table(nn_dx, file=paste0(newOut, "/knn_", newID, "_", Sys.Date(), ".txt"), row.names=T, quote=F, sep='\t')

    ##-----------------------------
    ## archetypal analysis: predict archetypes on new sample
    ##-----------------------------

    pred_aa <- data.frame(stats::predict(aa_model, new_scores[,colnames(aa_model$archetypes)]))
    rownames(pred_aa) <- newID
    pred_aa <- round(pred_aa * 100, 3)
    
    new_scores$aa_cluster <- gsub("X", "", colnames(pred_aa)[max.col(pred_aa)])
    new_scores$Dx <- "new"
    new_scores$Dx_simple <- "new"
    
    mscores_aa_all <- rbind(dx_ref, new_scores[,colnames(dx_ref)])

    pca.aa <- merge(pca.df, mscores_aa_all[,c("ID", "aa_cluster")], by="ID")
    pca.aa$aa_cluster <- as.factor(pca.aa$aa_cluster)

    pca.aa[pca.aa$ref=="new",]

    aa_cols=c("black", "firebrick", "lightgray", "blue3", "#009E73")

    pca_aa <- ggplot() +
        scale_fill_manual(values=aa_cols) +
        #geom_point(data=pca.aa, aes(Dim.1, Dim.2, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) +
        geom_point(data=pca.aa[pca.aa$ref=="ref",], aes(Dim.1, Dim.2, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) +
        geom_point(data=pca.aa[pca.aa$ref=="new",], aes(Dim.1, Dim.2), shape=23, size=4, alpha=1, fill="orange") +
        #geom_text_repel(data=pca.aa[pca.aa$ref=="new",], aes(Dim.1, Dim.2, label=ID), size=4, colour="orange") +
        xlab(paste("PC1 ", round(resPCA$eig[1,2]),"% of variance",sep="")) +
        ylab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
        theme(legend.position="top", legend.title=element_blank(),
              panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="white"),
              axis.text=element_text(size=12),
              axis.text.x=element_text(size=12, angle=0),
              axis.title.x=element_text(face="bold", size=12, angle=0),
              axis.title.y=element_text(face="bold", size=12),
              panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))

    pca_aa_2_3 <- ggplot() +
        scale_fill_manual(values=aa_cols) +
        geom_point(data=pca.aa[pca.aa$ref=="ref",], aes(Dim.2, Dim.3, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) +
        geom_point(data=pca.aa[pca.aa$ref=="new",], aes(Dim.2, Dim.3), shape=23, size=3, alpha=1, fill="orange") +
        #geom_text_repel(data=pca.aa[pca.aa$ref=="new",], aes(Dim.2, Dim.3, label=ID), size=3, colour="orange") +
        xlab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
        ylab(paste("PC3 ", round(resPCA$eig[3,2]),"% of variance",sep="")) +
        theme(legend.position="top", legend.title=element_blank(),
              panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="white"),
              axis.text=element_text(size=12),
              axis.text.x=element_text(size=12, angle=0),
              axis.title.x=element_text(face="bold", size=12, angle=0),
              axis.title.y=element_text(face="bold", size=12),
              panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    ##---------------
    ## Dx by cluster
    ## banff score by cluster
    # aa.df <- pca.aa[,c("ID", "aa_cluster", "Dx")]
    # #aa.df <- aa.df[aa.df$Dx!="new",]
    #
    # aa.1.df <- droplevels(aa.df[aa.df$aa_cluster=="1",])
    # aa.2.df <- droplevels(aa.df[aa.df$aa_cluster=="2",])
    # aa.3.df <- droplevels(aa.df[aa.df$aa_cluster=="3",])
    # aa.4.df <- droplevels(aa.df[aa.df$aa_cluster=="4",])
    #
    # ## Dx by cluster
    # k1.df <- data.frame(table(aa.1.df$Dx))
    # k1.df$k <- "1"
    # colnames(k1.df) <- c("Dx", "total", "cluster")
    #
    # k2.df <- data.frame(k2=table(aa.2.df$Dx))
    # k2.df$k <- "2"
    # colnames(k2.df) <- c("Dx", "total", "cluster")
    #
    # k3.df <- data.frame(k3=table(aa.3.df$Dx))
    # k3.df$k <- "3"
    # colnames(k3.df) <- c("Dx", "total", "cluster")
    #
    # k4.df <- data.frame(k4=table(aa.4.df$Dx))
    # k4.df$k <- "4"
    # colnames(k4.df) <- c("Dx", "total", "cluster")
    #
    # ## Dx by cluster
    # cluster_table <- dplyr::bind_rows(k1.df, k2.df, k3.df, k4.df)
    #write.table(cluster_table, file=paste0(newOut, "/aa_cluster_", newID, "_", Sys.Date(), ".txt"), row.names=F, quote=F, sep='\t')

    ##-------------------------------------------
    ## output files and plots for markdown report
    return(list(new_scores=new_scores_pct,
              ref_scores=ref.score.quantiles,
              knn_dx=nn_dx,
              knn_banff=nn_banff,
              pca_1_2=pca_new_1_2,
              pca_2_3=pca_new_2_3,
              aa_cluster_new=pred_aa, #new sample cluster probs
              pca_archetype=pca_aa,
              #aa_cluster_table=cluster_table,
              pathways=pathway_table,
              pathway_radar=pathway_radar,
              cell_types=cell_type_table,
              cell_type_radar=cell_type_radar,
              bkv_plot=bkv_boxplot,
              bkv_stats=bkv_tab)
          )

}
