
## predict probability of rejection on new biopsy
source('../scripts/BHOT.R')

##-----------------------------------------------------
## load models
##-----------------------------------------------------
modelPath='../model_data/BHOT/models/'

aa_model <- get(load(paste0(modelPath, 'aa_model.rda')))

## diagnosis
normal.model <- get(load(paste0(modelPath, 'normal_model.RData')))
amr.model <- get(load(paste0(modelPath, 'amr_model.RData')))
tcmr.model <- get(load(paste0(modelPath, 'tcmr_model.RData')))
atcmr.model <- get(load(paste0(modelPath, 'atcmr_model.RData')))
catcmr.model <- get(load(paste0(modelPath, 'catcmr_model.RData')))
ifta.model <- get(load(paste0(modelPath, 'ifta_model.RData')))

## banff lesions
## ordinal
g_score.model <- get(load(paste0(modelPath, 'g_ordinal_model.RData')))
ptc_score.model <- get(load(paste0(modelPath, 'ptc_ordinal_model.RData')))
i_score.model <- get(load(paste0(modelPath, 'i_ordinal_model.RData')))
t_score.model <- get(load(paste0(modelPath, 't_ordinal_model.RData')))
## binary
g0_score.model <- get(load(paste0(modelPath, 'g0_score_model.RData')))
ptc0_score.model <- get(load(paste0(modelPath, 'ptc0_score_model.RData')))
i1_score.model <- get(load(paste0(modelPath, 'i1_score_model.RData')))
t1_score.model <- get(load(paste0(modelPath, 't1_score_model.RData')))

cg0_score.model <- get(load(paste0(modelPath, 'cg0_score_model.RData')))
v0_score.model <- get(load(paste0(modelPath, 'v0_score_model.RData')))
iifta0_score.model <- get(load(paste0(modelPath, 'iifta0_score_model.RData')))
ci1_score.model <- get(load(paste0(modelPath, 'ci1_score_model.RData')))
ct1_score.model <- get(load(paste0(modelPath, 'ct1_score_model.RData')))
cv1_score.model <- get(load(paste0(modelPath, 'cv1_score_model.RData')))

## import refSet molecular scores and Dx
mscores_all_file <- '../model_data/BHOT/counts/refset_molecular_scores_all.txt'
mscores_aa_file <- '../model_data/BHOT/counts/refset_molecular_scores_aa.txt'

mscores_ref <- read.table(mscores_all_file, header=TRUE, sep='\t')
rownames(mscores_ref) <- mscores_ref$ID
colnames(mscores_ref)[colnames(mscores_ref)=="amr"]<-"AMR"
colnames(mscores_ref)[colnames(mscores_ref)=="tcmr"]<-"TCMR"
colnames(mscores_ref)[colnames(mscores_ref)=="atcmr"]<-"aTCMR"
colnames(mscores_ref)[colnames(mscores_ref)=="catcmr"]<-"caTCMR"
colnames(mscores_ref)[colnames(mscores_ref)=="ati"]<-"ATI"
colnames(mscores_ref)[colnames(mscores_ref)=="ifta"]<-"IFTA"

## import Banff scores
banff_ref <- read.table('../model_data/BHOT/counts/refset_banff_scores_all.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)

## Define histology diagnosis categories
dx_normal <- c("Normal or minimal changes", "Pristine", "No specific diagnosis")
dx_amr <- c("Active AMR", "Chronic (+/- active) AMR")
dx_tcmr <- c("Acute TCMR", "Chronic active TCMR")
dx_ati <- c("Acute tubular injury")
dx_ifta <- c("Isolated IFTA") #ie. IFTA w/ no rejection
dx_no_rejection <- unique(mscores_ref$Dx)[!unique(mscores_ref$Dx) %in% c(dx_amr, dx_tcmr)]

##-----------------------------------------------------
## import BHOT gene annotations
##-----------------------------------------------------
bhot_annot <- read.csv('../static/BHOT_annotations_il6.csv', check.names=FALSE, header=TRUE)
## exclude control and viral genes (n=12+4)
rm_genes <- bhot_annot[bhot_annot$`Internal Reference Gene`=="+" | bhot_annot$`Viral Detection`=="+","Gene"]
bhot_annot <- bhot_annot[!bhot_annot$Gene %in% rm_genes,]

endats <- read.table('../static/ENDAT_genes.txt', check.names=FALSE, header=FALSE)
endats <- endats[endats$V1 %in% bhot_annot$Gene,]

bhot_cell_types <- c("B-cells", "Macrophages", "T-cells", "NK cells")

bhot_pathways <- c("B-cell Receptor Signaling", "Chemokine Signaling", "Complement System", "IL6 Signaling", "MAPK", "mTOR", "NF-kappaB Signaling",
                   "Th1 Differentiation", "Th17 Differentiation", "Th2 Differentiation", "TNF Family Signaling", "Treg Differentiation",
                   "Type I Interferon Signaling", "Type II Interferon Signaling")

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

## manually exclude pediatric biospy RCCs for now
pedRCC <- read.table('../static/ped_rcc_ids.txt', header=FALSE)
refRCC <- refRCC[!basename(refRCC) %in% pedRCC[,1]]

#newRCC='../test_files/test.RCC'
#outPath='~/Downloads/'
#BHOTpred(newRCC, outPath)

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
      cat(">>Outputting results for sample:", newID, "\n")
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
    ## rename newID if also exists in refset sample names: rownames(mscores_ref)
    #samp_ind <- grep(newID, colnames(countTable))
    #if (length(samp_ind) > 1) {
    #  colnames(countTable)[samp_ind][1]<-newID
    #  colnames(countTable)[samp_ind][2]<-"histomx-new"
    #}
    ## keep only samples used in classifiers: eventually exclude unused RCC files from final repo
    #countTable <- countTable[,colnames(countTable) %in% c("CodeClass", "Name", "Accession", rownames(mscores_ref), newID)]
    
    ##------------------
    ## output raw counts for all reference samples used in classifiers
    #write.table(countTable, file='../static/refset_raw_counts.txt', quote=F, row.names=F, sep='\t')
    
    ## load raw count table, parse newRCC, and merge counts
    ns.data <- read.table('../model_data/BHOT/counts/refset_raw_counts.txt', sep='\t', header=TRUE, check.names=FALSE)
    ns.new <- parseRCC(newRCC)
    countTable <- merge(ns.data, ns.new$counts, by=c("CodeClass", "Name", "Accession"))
    rownames(countTable) <- countTable$Name
    
    ## TODO: ruvseq
    ##------------------------
    ## Normalization: no background correction
    x <- countTable[,!colnames(countTable) %in% c("CodeClass", "Name", "Accession")]

    ## 01 arithmetic mean of the geometric means of positive controls
    posTab <- t(x[rownames(x) %in% countTable[countTable$CodeClass=="Positive","Name"],])
    pos.sample <- apply(posTab, MARGIN=1, FUN=geoMean);
    pos.norm.factors <- data.frame(ID=names(pos.sample), mean=mean(pos.sample) / pos.sample)
    pos.norm.factors$mean <- round(pos.norm.factors$mean, 2)
    rownames(pos.norm.factors) <- NULL
    
    ## multiply normalization factor by raw counts
    x <- t(apply(x, MARGIN = 1, FUN = '*', pos.norm.factors$mean));

    ## 02 SampleContent (normalize to HK genes to account for sample or RNA content ie. pipetting fluctuations)
    ## Normalize by substracting geometric mean of housekeeping genes from each endogenous gene
    
    hk_genes=countTable[countTable$CodeClass=="Housekeeping","Name"]
    
    ## calculate the normalization factor: arithmetic mean of the geometric means of positive controls
    rna.content <- apply(x[rownames(x) %in% hk_genes,], MARGIN=2, FUN=geoMean);
    hk.norm.factor <- mean(rna.content) / rna.content
    
    x.norm <- t(apply(x, MARGIN = 1, FUN = '*', hk.norm.factor));
    x.norm <- log2(x.norm + 1);
    
    ## return endogenous probes
    ns.norm <- x.norm[rownames(x.norm) %in% countTable[countTable$CodeClass=="Endogenous","Name"],]
    ##------------------------
    
    ## new sample; proper delimiter for ordinal regression
    new.ns.norm.or <- data.frame(counts=ns.norm[,newID])
    colnames(new.ns.norm.or) <- newID
    
    ## TODO: prevent gene delimiter modification in model output
    rownames(ns.norm) <- gsub("-", ".", rownames(ns.norm))
    rownames(ns.norm) <- gsub("/", ".", rownames(ns.norm))
    rownames(ns.norm) <- gsub(" ", ".", rownames(ns.norm))
  
    ## new sample(s) normalized with refSet
    new.ns.norm <- data.frame(counts=ns.norm[,newID])
    colnames(new.ns.norm) <- newID
    
    ##--------------------------
    ## BKV expression: new biopsy v. normal and BKV refset samples
    ##--------------------------
    bk_genes <- c("BK..large.T.Ag", "BK..VP1")
    norm_bx_ids <- rownames(mscores_ref[mscores_ref$Dx %in% dx_normal,])
    bk_bx_ids <- rownames(mscores_ref[mscores_ref$Dx %in% c("NBKv"),])
      
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
    dat <- data.frame(t(bk_counts), check.names=TRUE)
    dat$group <- ifelse(rownames(dat) %in% norm_bx_ids, "normal",
                        ifelse(rownames(dat) %in% bk_bx_ids, "BKV",
                               ifelse(rownames(dat) %in% newID, newID, rownames(dat))))
    dat$group <- factor(dat$group, levels=c("normal", "BKV", newID))
    suppressWarnings({dat.m <- reshape2::melt(dat, id.vars="group")})
    bkv_boxplot <- ggplot(dat.m, aes(x = forcats::fct_rev(group), y = value, fill=group)) + geom_boxplot() +
      xlab("") + ylab("normalized expression") +
      scale_fill_manual(values=c("blue3",  "lightsteelblue", "orange")) + 
      theme(legend.position="none", panel.border=element_rect(colour="gray", fill=NA, size=1),
            axis.text=element_text(size=16, color="black"), axis.title=element_text(size=16, color="black"),
            panel.background=element_blank(), panel.grid.minor=element_line(colour="gray"))
    
    ## Wilcoxon rank sum test
    ## BK VP1
    vp1_new_norm_pval <- wilcox.test(dat[dat$group %in% c(newID),"BKVP1"], dat[dat$group %in% c("normal"),"BKVP1"], paired=FALSE, alternative="greater")$p.value
    vp1_new_bkv_pval <- wilcox.test(dat[dat$group %in% c(newID),"BKVP1"], dat[dat$group %in% c("BKV"),"BKVP1"], paired=FALSE, alternative="less")$p.value
    vp1_bkv_norm_pval <- wilcox.test(dat[dat$group %in% c("BKV"),"BKVP1"], dat[dat$group %in% c("normal"),"BKVP1"], paired=FALSE, alternative="greater")$p.value
    # BK large T Ag
    tag_new_norm_pval <- wilcox.test(dat[dat$group %in% c(newID),"BKlargeTAg"], dat[dat$group %in% c("normal"),"BKlargeTAg"], paired=FALSE, alternative="greater")$p.value
    tag_new_bkv_pval <- wilcox.test(dat[dat$group %in% c(newID),"BKlargeTAg"], dat[dat$group %in% c("BKV"),"BKlargeTAg"], paired=FALSE, alternative="less")$p.value
    tag_bkv_norm_pval <- wilcox.test(dat[dat$group %in% c("BKV"),"BKlargeTAg"], dat[dat$group %in% c("normal"),"BKlargeTAg"], paired=FALSE, alternative="greater")$p.value
    
    bkv_tab$new_normal_pval <- formatC(c(vp1_new_norm_pval, tag_new_norm_pval), format="e", digits=3)
    bkv_tab$new_bkv_pval <- formatC(c(vp1_new_bkv_pval, tag_new_bkv_pval), format="e", digits=3)
    bkv_tab$bkv_normal_pval <- formatC(c(vp1_bkv_norm_pval, tag_bkv_norm_pval), format="e", digits=3)
    
    ##--------------------------
    ## signaling pathways
    ##--------------------------
    norm_bx_ids <- rownames(mscores_ref[mscores_ref$Dx %in% dx_normal,])
    #norm_bx_ids <- rownames(mscores_ref[mscores_ref$Dx %in% dx_no_rejection,])
  
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
    ## predict Dx/Banff score probabilities
    ##--------------------------
    ## predict normal prob for new sample(s)
    #normal.genes <- rownames(normal.model$scaling) #LDA
    normal.genes <- names(normal.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% normal.genes]
  
    #new.normal.pred <- predict(normal.model, newdata=new.ge, type="prob")$posterior #LDA
    new.normal.pred <- predict(normal.model, newdata=data.frame(t(new.ge)), type="response")
  
    new.normal.pred <- data.frame(new.normal.pred)
    #new.normal.pred$ID <- rownames(new.normal.pred)
    new.normal.pred$ID <- newID
    #new.normal.pred$other <- NULL #LDA
    colnames(new.normal.pred)<-c("normal", "ID")
  
    ##--------------------------
    ## predict AMR prob for new sample(s)
    #amr.genes <- rownames(amr.model$scaling) #LDA
    amr.genes <- names(amr.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% amr.genes]
  
    #new.amr.pred <- predict(amr.model, newdata=new.ge, type="prob")$posterior #LDA
    new.amr.pred <- predict(amr.model, newdata=data.frame(t(new.ge)), type="response")
  
    new.amr.pred <- data.frame(new.amr.pred)
    #new.amr.pred$ID <- rownames(new.amr.pred)
    new.amr.pred$ID <- newID
    #new.amr.pred$other <- NULL #LDA
    colnames(new.amr.pred)<-c("amr", "ID")
  
    ##--------------------------
    ## predict TCMR prob for new sample(s)
    #tcmr.genes <- rownames(tcmr.model$scaling) #LDA
    tcmr.genes <- names(tcmr.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% tcmr.genes]
  
    #new.tcmr.pred <- predict(tcmr.model, newdata=new.ge, type="prob")$posterior #LDA
    new.tcmr.pred <- predict(tcmr.model, newdata=data.frame(t(new.ge)), type="response")
  
    new.tcmr.pred <- data.frame(new.tcmr.pred)
    #new.tcmr.pred$ID <- rownames(new.tcmr.pred)
    new.tcmr.pred$ID <- newID
    #new.tcmr.pred$other <- NULL #LDA
    colnames(new.tcmr.pred)<-c("tcmr", "ID")
    
    ##--------------------------
    ## predict acute TCMR prob for new sample(s)
    #atcmr.genes <- rownames(atcmr.model$scaling) #LDA
    atcmr.genes <- names(atcmr.model$coefficients)[-1] #LR
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% atcmr.genes]
    
    #new.atcmr.pred <- predict(atcmr.model, newdata=new.ge, type="prob")$posterior #LDA
    new.atcmr.pred <- predict(atcmr.model, newdata=data.frame(t(new.ge)), type="response")
    
    new.atcmr.pred <- data.frame(new.atcmr.pred)
    #new.atcmr.pred$ID <- rownames(new.atcmr.pred)
    new.atcmr.pred$ID <- newID
    #new.atcmr.pred$other <- NULL #LDA
    colnames(new.atcmr.pred)<-c("atcmr", "ID")
  
    ##--------------------------
    ## predict acute TCMR prob for new sample(s)
    #catcmr.genes <- rownames(catcmr.model$scaling) #LDA
    catcmr.genes <- names(catcmr.model$coefficients)[-1] #LR
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% catcmr.genes]
    
    #new.catcmr.pred <- predict(catcmr.model, newdata=new.ge, type="prob")$posterior #LDA
    new.catcmr.pred <- predict(catcmr.model, newdata=data.frame(t(new.ge)), type="response")
    
    new.catcmr.pred <- data.frame(new.catcmr.pred)
    #new.catcmr.pred$ID <- rownames(new.catcmr.pred)
    new.catcmr.pred$ID <- newID
    #new.catcmr.pred$other <- NULL #LDA
    colnames(new.catcmr.pred)<-c("catcmr", "ID")
  
    ##--------------------------
    ## predict IFTA prob for new sample(s)
    #ifta.genes <- rownames(ifta.model$scaling) #LDA
    ifta.genes <- names(ifta.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% ifta.genes]
  
    #new.ifta.pred <- predict(ifta.model, newdata=new.ge, type="prob")$posterior #LDA
    new.ifta.pred <- predict(ifta.model, newdata=data.frame(t(new.ge)), type="response")
  
    new.ifta.pred <- data.frame(new.ifta.pred)
    #new.ifta.pred$ID <- rownames(new.ifta.pred)
    new.ifta.pred$ID <- newID
    #new.ifta.pred$other <- NULL #LDA
    colnames(new.ifta.pred)<-c("ifta", "ID")
  
    ##--------------------------
    ## predict prob g>0 for new sample(s)
    #g0.genes <- rownames(g0_score.model$scaling) #LDA
    g0.genes <- names(g0_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% g0.genes]
  
    #new.pred <- predict(g0_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(g0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    #new.pred$ID <- rownames(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.g0.pred <- new.pred
    colnames(new.g0.pred)<-c("g0_score", "ID")
  
    ##--------------------------
    ## predict prob ptc>0 for new sample(s)
    #ptc0.genes <- rownames(ptc0_score.model$scaling) #LDA
    ptc0.genes <- names(ptc0_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% ptc0.genes]
  
    #new.pred <- predict(ptc0_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(ptc0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    #new.pred$ID <- rownames(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ptc0.pred <- new.pred
    colnames(new.ptc0.pred) <- c("ptc0_score", "ID")
  
    ##--------------------------
    ## predict prob i>1 for new samples
    #i1.genes <- rownames(i1_score.model$scaling) #LDA
    i1.genes <- names(i1_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% i1.genes]
  
    #new.pred <- predict(i1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(i1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    #new.pred$ID <- rownames(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.i1.pred <- new.pred
    colnames(new.i1.pred) <- c("i1_score", "ID")
  
    ##--------------------------
    ## predict prob t>1 for new samples
    #t1.genes <- rownames(t1_score.model$scaling) #LDA
    t1.genes <- names(t1_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% t1.genes]
  
    #new.pred <- predict(t1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(t1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    #new.pred$ID <- rownames(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.t1.pred <- new.pred
    colnames(new.t1.pred) <- c("t1_score", "ID")
  
    ##--------------------------
    ## ordinal regression
    
    ## g score
    g.genes <- names(g_score.model$coefficients)[-c(1:3)]
    g.genes <- gsub("`", '', g.genes)
    
    new.ge <- new.ns.norm.or
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% g.genes]
    
    new.g.pred <- data.frame(g=predict(g_score.model, newdata=new.ge, type="prob"))
    colnames(new.g.pred) <- gsub(".fit.", "", colnames(new.g.pred))
    new.g.pred$ID <- newID
    
    ## ptc score
    ptc.genes <- names(ptc_score.model$coefficients)[-c(1:3)]
    ptc.genes <- gsub("`", '', ptc.genes)
    
    new.ge <- new.ns.norm.or
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% ptc.genes]
    
    new.ptc.pred <- data.frame(ptc=predict(ptc_score.model, newdata=new.ge, type="prob"))
    colnames(new.ptc.pred) <- gsub(".fit.", "", colnames(new.ptc.pred))
    new.ptc.pred$ID <- newID
    
    ## i score
    i.genes <- names(i_score.model$coefficients)[-c(1:3)]
    i.genes <- gsub("`", '', i.genes)
    
    new.ge <- new.ns.norm.or
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% i.genes]
    
    new.i.pred <- data.frame(i=predict(i_score.model, newdata=new.ge, type="prob"))
    colnames(new.i.pred) <- gsub(".fit.", "", colnames(new.i.pred))
    new.i.pred$ID <- newID
    
    ## t score
    t.genes <- names(t_score.model$coefficients)[-c(1:3)]
    t.genes <- gsub("`", '', t.genes)
    
    new.ge <- new.ns.norm.or
    new.ge <- data.frame(t(new.ge), check.names=FALSE)
    new.ge <- new.ge[,colnames(new.ge) %in% t.genes]
    
    new.t.pred <- data.frame(t=predict(t_score.model, newdata=new.ge, type="prob"))
    colnames(new.t.pred) <- gsub(".fit.", "", colnames(new.t.pred))
    new.t.pred$ID <- newID
    
    ##--------------------------
    ## g, ptc, i, t binary classifiers
    
    ##---------------
    #g0.genes <- rownames(g0_score.model$scaling) #LDA
    g0.genes <- names(g0_score.model$coefficients)[-1] #LR 
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% g0.genes]
    
    #new.pred <- predict(g0_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(g0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
    
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.g0.pred <- new.pred
    colnames(new.g0.pred) <- c("g0_score", "ID")
    
    ##---------------
    #ptc0.genes <- rownames(ptc0_score.model$scaling) #LDA
    ptc0.genes <- names(ptc0_score.model$coefficients)[-1] #LR 
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% ptc0.genes]
    
    #new.pred <- predict(ptc0_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(ptc0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
    
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ptc0.pred <- new.pred
    colnames(new.ptc0.pred) <- c("ptc0_score", "ID")
    
    ##---------------
    #i1.genes <- rownames(i1_score.model$scaling) #LDA
    i1.genes <- names(i1_score.model$coefficients)[-1] #LR 
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% i1.genes]
    
    #new.pred <- predict(i1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(i1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
    
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.i1.pred <- new.pred
    colnames(new.i1.pred) <- c("i1_score", "ID")
    
    ##---------------
    #t1.genes <- rownames(t1_score.model$scaling) #LDA
    t1.genes <- names(t1_score.model$coefficients)[-1] #LR 
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% t1.genes]
    
    #new.pred <- predict(t1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(t1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
    
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.t1.pred <- new.pred
    colnames(new.t1.pred) <- c("t1_score", "ID")
    
    ##--------------------------
    ## predict prob ci>1 for new sample(s)
    #ci1.genes <- rownames(ci1_score.model$scaling) #LDA
    ci1.genes <- names(ci1_score.model$coefficients)[-1] #LR 
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% ci1.genes]
  
    #new.pred <- predict(ci1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(ci1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ci1.pred <- new.pred
    colnames(new.ci1.pred) <- c("ci1_score", "ID")
  
    ##--------------------------
    ## predict prob ct>1 for new sample(s)
    #ct1.genes <- rownames(ct1_score.model$scaling) #LDA
    ct1.genes <- names(ct1_score.model$coefficients)[-1] #LR 
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% ct1.genes]
  
    #new.pred <- predict(ct1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(ct1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.ct1.pred <- new.pred
    colnames(new.ct1.pred) <- c("ct1_score", "ID")
  
    ##--------------------------
    ## predict prob cv>1 for new sample(s)
    #cv1.genes <- rownames(cv1_score.model$scaling) #LDA
    cv1.genes <- names(cv1_score.model$coefficients)[-1] #LR 
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% cv1.genes]
  
    #new.pred <- predict(cv1_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(cv1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.cv1.pred <- new.pred
    colnames(new.cv1.pred) <- c("cv1_score", "ID")
  
    ##--------------------------
    ## predict prob v>0 for new sample(s)
    #v0.genes <- rownames(v0_score.model$scaling) #LDA
    v0.genes <- names(v0_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% v0.genes]
  
    #new.pred <- predict(v0_score.model, newdata=new.ge, type="prob")$posterior #LDA
    new.pred <- predict(v0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.v0.pred <- new.pred
    colnames(new.v0.pred) <- c("v0_score", "ID")
  
    ##--------------------------
    ## predict prob cg>0 for new samples
    #cg0.genes <- rownames(cg0_score.model$scaling) #LDA
    cg0.genes <- names(cg0_score.model$coefficients)[-1] #LR
  
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% cg0.genes]

    #new.pred <- predict(cg0_score.model, newdata=new.ge, type="prob")$posterior
    new.pred <- predict(cg0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.cg0.pred <- new.pred
    colnames(new.cg0.pred) <- c("cg0_score", "ID")
  
    ##--------------------------
    ## predict prob i_IFTA>0 for new samples
    #iifta0.genes <- rownames(cg0_score.model$scaling) #LDA
    iifta0.genes <- names(iifta0_score.model$coefficients)[-1] #LR
    
    new.ge <- new.ns.norm
    new.ge <- t(new.ge)
    new.ge <- new.ge[,colnames(new.ge) %in% iifta0.genes]
    
    #new.pred <- predict(iifta0_score.model, newdata=new.ge, type="prob")$posterior
    new.pred <- predict(iifta0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
    
    new.pred <- data.frame(new.pred)
    new.pred$ID <- newID
    new.pred$low <- NULL
    new.iifta0.pred <- new.pred
    colnames(new.iifta0.pred) <- c("iifta0_score", "ID")
    
    ##-------------------------------------
    ## score table
    options(scipen = 999)
    join_list <- list(new.normal.pred, new.amr.pred, new.ifta.pred,
                      new.tcmr.pred, new.atcmr.pred, new.catcmr.pred, 
                      new.g.pred, new.ptc.pred, new.i.pred, new.t.pred, #ordinal
                      new.g0.pred, new.ptc0.pred, new.i1.pred, new.t1.pred, #binary
                      new.cg0.pred, new.iifta0.pred, new.v0.pred,
                      new.cv1.pred, new.ci1.pred, new.ct1.pred)
    tab <- Reduce(function(...) merge(..., all=TRUE), join_list)
    rownames(tab) <- tab$ID
  
    colnames(tab)[colnames(tab)=="amr"]<-"AMR"
    colnames(tab)[colnames(tab)=="tcmr"]<-"TCMR"
    colnames(tab)[colnames(tab)=="atcmr"]<-"aTCMR"
    colnames(tab)[colnames(tab)=="catcmr"]<-"caTCMR"
    colnames(tab)[colnames(tab)=="ifta"]<-"IFTA"
  
    mscores_new <- tab ## save scores before calculating %
    #tab <- data.frame(apply(tab[,-1], 2, function(x) {round(x, 3)*100})) #multiple samples
    tab <- data.frame(round(tab[-1], 3) *100) #single sample
    tab$ID <- rownames(tab)
    mscores_new_pct <- tab %>% dplyr::select(ID, everything())
  
    ## output score table
    if(saveFiles=="TRUE") {
      write.table(mscores_new_pct, file=paste0(newOut, "/molecular_score_table_", newID, "_", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
    }
  
    options(scipen = 0)

    ##-----------------------
    ## output IQR and median scores by Dx for reference biopsies
    scores_keep <- c("AMR", "TCMR", "aTCMR", "caTCMR", "IFTA", "normal",
                     "g0", "g1", "g2", "g3", "ptc0", "ptc1", "ptc2", "ptc3",
                     "i1", "i2", "i3", "t0", "t1", "t2", "t3",
                     "g0_score", "ptc0_score", "i1_score", "t1_score",
                     "cg0_score", "iifta0_score", "v0_score",
                     "ci1_score", "ct1_score", "cv1_score")
  
    ref.tab <- mscores_ref ## complete Dx + molecular scores
    #ref.tab <- ref.tab[apply(ref.tab[,-1], 1, function(x) { all(!is.na(x)) }),] #exclude samples without all molecular scores
  
    ## simplify Dx
    ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_amr, "AMR", ref.tab$Dx)
    ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_tcmr, "TCMR", ref.tab$Dx)
    ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_normal, "No specific Dx", ref.tab$Dx)
    ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_ifta, "IFTA", ref.tab$Dx)
  
    ## IQR = Q3 - Q1
    ## Q0= min and Q5=max
    ref.amr.quantiles <- apply(ref.tab[ref.tab$Dx=="AMR",scores_keep], 2, quantile, na.rm=TRUE)
    ref.tcmr.quantiles <- apply(ref.tab[ref.tab$Dx=="TCMR",scores_keep], 2, quantile, na.rm=TRUE)
    #ref.tcmr.quantiles <- apply(ref.tab[ref.tab$Dx %in% dx_tcmr,scores_keep], 2, quantile, na.rm=TRUE)
    ref.ifta.quantiles <- apply(ref.tab[ref.tab$Dx=="IFTA",scores_keep], 2, quantile, na.rm=TRUE)
    ref.normal.quantiles <- apply(ref.tab[ref.tab$Dx=="No specific Dx",scores_keep], 2, quantile, na.rm=TRUE)
  
    ref.amr.quantiles <- data.frame(apply(ref.amr.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.amr.quantiles$Dx <- "AMR"
    ref.amr.quantiles$Q <- rownames(ref.amr.quantiles)
  
    ref.tcmr.quantiles <- data.frame(apply(ref.tcmr.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.tcmr.quantiles$Dx <- "TCMR"
    ref.tcmr.quantiles$Q <- rownames(ref.tcmr.quantiles)
  
    ref.ifta.quantiles <- data.frame(apply(ref.ifta.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.ifta.quantiles$Dx <- "IFTA"
    ref.ifta.quantiles$Q <- rownames(ref.ifta.quantiles)
  
    ref.normal.quantiles <- data.frame(apply(ref.normal.quantiles, 2, function(x) {round(x, 3)*100}))
    ref.normal.quantiles$Dx <- "No specific Dx"
    ref.normal.quantiles$Q <- rownames(ref.normal.quantiles)
  
    ref.score.quantiles <- rbind(ref.amr.quantiles, ref.tcmr.quantiles, ref.ifta.quantiles, ref.normal.quantiles)
    #write.table(ref.score.quantiles, file=paste0(newOut, "/refset_score_quantiles", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
  
    ## resfset median scores (needed to reorder boxplots)
    ref.amr.median <- apply(ref.tab[ref.tab$Dx=="AMR",scores_keep], 2, median, na.rm=TRUE)
    ref.tcmr.median <- apply(ref.tab[ref.tab$Dx=="TCMR",scores_keep], 2, median, na.rm=TRUE)
    ref.ifta.median <- apply(ref.tab[ref.tab$Dx=="IFTA",scores_keep], 2, median, na.rm=TRUE)
    ref.normal.median <- apply(ref.tab[ref.tab$Dx=="No specific Dx",scores_keep], 2, median, na.rm=TRUE)
    #ref.scores <- rbind(ref.amr.median, ref.tcmr.median, ref.normal.median)
    #ref.scores <- data.frame(apply(ref.scores, 1, function(x) {round(x, 3)*100}))
    #colnames(ref.scores) <- gsub("ref.", "", colnames(ref.scores))
    #ref.scores$score <- rownames(ref.scores)
    ##write.table(ref.scores, file=paste0(newOut, "/refset_median_scores_", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
  
    ##----------------------------------------------
    ## boxplots of reference biopsy molecular scores
    ## AMR
    dat <- ref.tab[ref.tab$Dx=="AMR",]
    dat <- dat[,scores_keep]
    amr_order <- names(ref.amr.median[order(ref.amr.median, decreasing=F)])
    ## order by median scores (reorder not working)
    dat <- dat[,amr_order]
    colnames(dat) <- gsub("_score", "", colnames(dat))
    suppressWarnings({
      dat.m <- reshape2::melt(dat)
    })
    dat.m$value <- dat.m$value * 100
  
    #ggplot(data=dat.m, aes(x=reorder(variable, value, median, order=TRUE), y=value)) + 
    boxplot_amr <- ggplot(data=dat.m, aes(x=variable, y=value)) + 
        ggtitle("AMR reference biopsies") +
        ylab("molecular score (%)") + xlab("") +
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
              panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    #ggsave(paste0(newOut, "/amr_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_amr, device="pdf", width=8, height=6)
  
    ## TCMR
    dat <- ref.tab[ref.tab$Dx=="TCMR",]
    dat <- dat[,scores_keep]
    tcmr_order <- names(ref.tcmr.median[order(ref.tcmr.median, decreasing=F)])
    ## order by median scores (reorder not working)
    dat <- dat[,tcmr_order]
    colnames(dat) <- gsub("_score", "", colnames(dat))
    suppressWarnings({
      dat.m <- data.table::melt(dat)
    })
    dat.m$value <- dat.m$value * 100
    boxplot_atcmr <- ggplot(data=dat.m, aes(x=variable, y=value)) + 
        ggtitle("TCMR reference biopsies") +
        ylab("molecular score (%)") + xlab("") +
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
              panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    #ggsave(paste0(newOut, "/atcmr_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_tcmr, device="pdf", width=8, height=6)

    ## Normal
    dat <- ref.tab[ref.tab$Dx=="No specific Dx",]
    dat <- dat[,scores_keep]
    normal_order <- names(ref.normal.median[order(ref.normal.median, decreasing=F)])
    ## order by median scores (reorder not working)
    dat <- dat[,normal_order]
    colnames(dat) <- gsub("_score", "", colnames(dat))
    suppressWarnings({
      dat.m <- data.table::melt(dat)
    })
    dat.m$value <- dat.m$value * 100
  
    #ggplot(data=dat.m, aes(x=reorder(variable, value, median, order=TRUE), y=value)) + 
    boxplot_normal <- ggplot(data=dat.m, aes(x=variable, y=value)) + 
        ggtitle("No specific Dx reference biopsies") +
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
              panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    #ggsave(paste0(newOut, "/normal_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_normal, device="pdf", width=8, height=6)
  
    ##--------------------------
    ## combine refset and new molecular scores
    ## scores to input to PCA
    #pca_scores <- c("normal", "AMR", "aTCMR", "IFTA")
    pca_scores <- c("normal", "AMR", "TCMR", "IFTA")
    #pca_scores <- c("g0_score", "ptc0_score", "cg0_score", "i1_score", "t1_score")
    
    ## keep only Bx with all molecular scores
    mscores_pca <- mscores_ref[apply(mscores_ref[,pca_scores], 1, function(x) { all(!is.na(x)) }),]
    
    mscores_pca <- rbind(mscores_pca[,pca_scores], mscores_new[,pca_scores])
      
    resPCA <- FactoMineR::PCA(mscores_pca[,!colnames(mscores_pca) %in% "ID"], scale.unit=FALSE, ncp=5, graph=FALSE)
  
    pca.df <- as.data.frame(resPCA$ind$coord)
    pca.df$ID <- rownames(pca.df)
  
    pca.df <- merge(pca.df, mscores_ref[,c("ID", "Dx")], by="ID", all.x=TRUE)
    pca.df$Dx <- ifelse(pca.df$ID %in% mscores_new$ID, "new", pca.df$Dx)
  
    ## higlight new biopsy
    pca.df <- plyr::mutate(pca.df, ref=ifelse(pca.df$Dx!="new", "ref", "new"))
  
    ## simple Dx
    pca.df$Dx <- ifelse(pca.df$Dx %in% c("Chronic (+/- active) AMR"), "Chronic active AMR", pca.df$Dx)
    #pca.df$Dx <- ifelse(pca.df$Dx %in% dx_amr, "AMR", pca.df$Dx)
    #pca.df$Dx <- ifelse(pca.df$Dx %in% dx_tcmr, "TCMR", pca.df$Dx)
    pca.df$Dx <- ifelse(pca.df$Dx %in% dx_normal, "No specific Dx", pca.df$Dx)
    pca.df$Dx <- ifelse(pca.df$Dx %in% dx_ifta, "IFTA", pca.df$Dx)
  
    #pca.df <- pca.df[pca.df$Dx %in% c("AMR", "TCMR", "ATI", "IFTA", "No specific Dx", "new"),]
    pca.df <- pca.df[pca.df$Dx %in% c("Active AMR",  "Chronic active AMR", "Acute TCMR", "Chronic active TCMR", "IFTA", "No specific Dx", "new"),]
  
    #my_cols=c("firebrick", "#009E73", "gray80", "black", "blue3") #steelblue;dodgerblue
    my_cols=c("darkred", "blue3", "tomato", "dodgerblue", "#009E73", "lightsteelblue", "black")
    paste(levels(factor(pca.df$Dx)), my_cols)
    #table(pca.df$Dx)
  
    pca_new_1_2 <- ggplot() +
        scale_fill_manual(values=my_cols) +
        #geom_point(data=pca.df, aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) + 
        geom_point(data=pca.df[pca.df$ref=="ref",], aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=3, alpha=0.8) + 
        geom_point(data=pca.df[pca.df$ref=="new",], aes(Dim.1, Dim.2), shape=23, size=4, alpha=1, fill="orange") +
        #geom_text_repel(data=pca.df[pca.df$ref=="new",], aes(Dim.1, Dim.2, label=ID), size=2, colour="orange") +
        xlab(paste("PC1 ", round(resPCA$eig[1,2]),"% of variance",sep="")) +
        ylab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
        theme(legend.position="right", legend.title=element_blank(),
              panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
              panel.background=element_blank(), 
              axis.line=element_line(colour="white"),
              axis.text=element_text(size=12), 
              axis.text.x=element_text(face="bold", size=12, angle=0),
              axis.title.x=element_text(face="bold", size=12, angle=0), 
              axis.title.y=element_text(face="bold", size=12),
              panel.border=element_rect(colour="white", fill=NA, size=5))
  
    pca_new_2_3 <- ggplot() +
      scale_fill_manual(values=my_cols) +
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
    #knn_tab <- rbind(mscores_ref[,pca_scores], mscores_new[,pca_scores])
    knn_tab <- mscores_pca
    knn_tab$ID <- rownames(knn_tab)
    knn_tab <- merge(knn_tab, pca.df[,c("ID", "Dx", "ref")], by="ID")
  
    train_ref <- knn_tab[knn_tab$ref %in% "ref",c("normal", "AMR", "TCMR", "IFTA", "ID", "Dx")]
    test_new <- knn_tab[knn_tab$ref %in% "new",c("normal", "AMR", "TCMR", "IFTA")]
    #train_ref <- knn_tab[knn_tab$ref=="ref",c("cg", "g", "ptc", "i", "t", "ID", "Dx")]
    #test_new <- knn_tab[knn_tab$ref=="new",c("cg0", "g0", "ptc0", "i1", "t1")]
  
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
    banff_df <- merge(nn_df, banff_ref[,c("ID", "g_score", "ptc_score", "cg_score", "i_score", "t_score")], by="ID")

    g_df <- data.frame(g=table(banff_df$g_score, useNA="no"))
    ptc_df <- data.frame(ptc=table(banff_df$ptc_score, useNA="no"))
    cg_df <- data.frame(cg=table(banff_df$cg_score, useNA="no"))
    i_df <- data.frame(i=table(banff_df$i_score, useNA="no"))
    t_df <- data.frame(t=table(banff_df$t_score, useNA="no"))
    
    score_df <- data.frame(score=c("0", "1", "2", "3"))
    score_df <- merge(score_df, g_df, by.x="score", by.y="g.Var1", all=TRUE)
    score_df <- merge(score_df, ptc_df, by.x="score", by.y="ptc.Var1", all=TRUE)
    score_df <- merge(score_df, cg_df, by.x="score", by.y="cg.Var1", all=TRUE)
    score_df <- merge(score_df, i_df, by.x="score", by.y="i.Var1", all=TRUE)
    score_df <- merge(score_df, t_df, by.x="score", by.y="t.Var1", all=TRUE)
    
    score_df[is.na(score_df)]<-0
    
    rownames(score_df) <- score_df$score
    nn_banff <- data.frame(round(apply(score_df[,-1], 2, function(x) {x/sum(x)*100}), 0), check.names=F) #% of Bx with given score by lesion
    colnames(nn_banff) <- gsub(".Freq", "", colnames(nn_banff))
    
    ##----------------
    ## mean molecular score of nearest neighbors
    ## TODO: if not all Dx in nn_dx table
    #nn_scores <- knn_tab[knn_tab$ref=="ref",c("ID", "normal", "AMR", "TCMR", "ATI")]
    #nn_scores <- nn_scores[nn_scores$ID %in% nn_df$ID,]
  
    #nn_means <- apply(nn_scores[,-1], 2, function(x) { mean(as.numeric(x)) })
    #nn_means <- round(nn_means * 100, 2) #convert to %
    #nn_means <- data.frame(t(nn_means))
    #rownames(nn_means) <- "mean_score_knn"
  
    #kable(nn_means, row.names=TRUE, align='l', caption="", format="html") %>% kableExtra::kable_styling(position="left")
    #nn_dx$mean_score_knn <- t(nn_means)
    #write.table(nn_dx, file=paste0(newOut, "/knn_", newID, "_", Sys.Date(), ".txt"), row.names=T, quote=F, sep='\t')
  
    ##-----------------------------
    ## archetypal analysis: predict archetypes on unseen molecular scores
    ##-----------------------------
    ## import cluster assignment for ref biopsies
    # mscores_aa <- read.table(mscores_aa_file, header=TRUE, sep='\t')
    # 
    # mscores_new <- dplyr::rename(mscores_new, amr=AMR)
    # mscores_new <- dplyr::rename(mscores_new, tcmr=TCMR)
    # mscores_new <- dplyr::rename(mscores_new, ifta=IFTA)
    # 
    # pred_aa <- data.frame(predict(aa_model, mscores_new[,colnames(aa_model$archetypes)]))
    # rownames(pred_aa) <- newID
    # mscores_new$aa_cluster <- gsub("X", "", colnames(pred_aa)[max.col(pred_aa)])
    # mscores_new$Dx <- "new"
    # 
    # mscores_aa_all <- rbind(mscores_aa, mscores_new[,colnames(mscores_aa)])
    # 
    # pca.aa <- merge(pca.df, mscores_aa_all[,c("ID", "aa_cluster")], by="ID")
    # pca.aa$aa_cluster <- as.factor(pca.aa$aa_cluster)
    # 
    # pca.aa[pca.aa$ref=="new",]
    # 
    # aa_cols=c("black", "firebrick", "gray80", "#009E73", "firebrick")
    # 
    # pca_aa <- ggplot() +
    #     scale_fill_manual(values=aa_cols) +
    #     #geom_point(data=pca.aa, aes(Dim.1, Dim.2, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) + 
    #     geom_point(data=pca.aa[pca.aa$ref=="ref",], aes(Dim.1, Dim.2, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) + 
    #     geom_point(data=pca.aa[pca.aa$ref=="new",], aes(Dim.1, Dim.2), shape=23, size=4, alpha=1, fill="orange") +
    #     #geom_text_repel(data=pca.aa[pca.aa$ref=="new",], aes(Dim.1, Dim.2, label=ID), size=4, colour="orange") +
    #     xlab(paste("PC1 ", round(resPCA$eig[1,2]),"% of variance",sep="")) +
    #     ylab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
    #     theme(legend.position="top", legend.title=element_blank(),
    #           panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
    #           panel.background=element_blank(), 
    #           axis.line=element_line(colour="white"),
    #           axis.text=element_text(size=12), 
    #           axis.text.x=element_text(size=12, angle=0),
    #           axis.title.x=element_text(face="bold", size=12, angle=0), 
    #           axis.title.y=element_text(face="bold", size=12),  
    #           panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    # 
    # pca_aa_2_3 <- ggplot() +
    #     scale_fill_manual(values=aa_cols) +
    #     geom_point(data=pca.aa[pca.aa$ref=="ref",], aes(Dim.2, Dim.3, fill=aa_cluster), shape=21, color="gray", size=4, alpha=0.7) + 
    #     geom_point(data=pca.aa[pca.aa$ref=="new",], aes(Dim.2, Dim.3), shape=23, size=3, alpha=1, fill="orange") +
    #     #geom_text_repel(data=pca.aa[pca.aa$ref=="new",], aes(Dim.2, Dim.3, label=ID), size=3, colour="orange") +
    #     xlab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
    #     ylab(paste("PC3 ", round(resPCA$eig[3,2]),"% of variance",sep="")) +
    #     theme(legend.position="top", legend.title=element_blank(),
    #           panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
    #           panel.background=element_blank(), 
    #           axis.line=element_line(colour="white"),
    #           axis.text=element_text(size=12), 
    #           axis.text.x=element_text(size=12, angle=0),
    #           axis.title.x=element_text(face="bold", size=12, angle=0), 
    #           axis.title.y=element_text(face="bold", size=12),  
    #           panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
    #if (saveFiles=="TRUE") {
      #ggsave(paste0(newOut, "/pca_archetypes_pc2_3_", newID, "_", Sys.Date(), ".pdf"), plot=last_plot(), device="pdf", width=7, height=6)
    #}
    
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
    # aa.5.df <- droplevels(aa.df[aa.df$aa_cluster=="5",])
    # aa.6.df <- droplevels(aa.df[aa.df$aa_cluster=="6",])
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
    # k5.df <- data.frame(k5=table(aa.5.df$Dx))
    # k5.df$k <- "5"
    # colnames(k5.df) <- c("Dx", "total", "cluster")
    # 
    # k6.df <- data.frame(k6=table(aa.6.df$Dx))
    # k6.df$k <- "6"
    # colnames(k6.df) <- c("Dx", "total", "cluster")
    # 
    # ## Dx by cluster
    # cluster_table <- dplyr::bind_rows(k1.df, k2.df, k3.df, k4.df, k5.df, k6.df)
    #write.table(cluster_table, file=paste0(newOut, "/aa_cluster_", newID, "_", Sys.Date(), ".txt"), row.names=F, quote=F, sep='\t')
  
    ##-------------------------------------------
    ## output files and plots for markdown report
    return(list(score_table=mscores_new_pct,
              ref_scores=ref.score.quantiles,
              knn_dx=nn_dx,
              knn_banff=nn_banff,
              pca_1_2=pca_new_1_2,
              pca_2_3=pca_new_2_3,
              #aa_cluster_table=cluster_table,
              #aa_cluster_new=pred_aa, #new sample cluster probs
              #pca_archetype=pca_aa,
              pathways=pathway_table,
              pathway_radar=pathway_radar,
              cell_types=cell_type_table,
              cell_type_radar=cell_type_radar,
              bkv_plot=bkv_boxplot,
              bkv_stats=bkv_tab)
          )
  
}

##--------------------------   
## barplots of banff scores by cluster
## k_df is a table with ID, aa_cluster, Dx + banff lesion scores
plot_k <- function(k_df, cluster, index) {
  
  df.tab <- mapply(table, k_df[,c("cg", "g", "ptc", "i", "t", "v", "ci", "ct", "cv", "ah")])
  
  if (isTRUE(index)) {
    ##index using largest vector as number of rows
    df.banff <- data.frame(sapply(df.tab, '[', seq(max(lengths(df.tab))))) 
  } else {
    df.banff <- df.tab
  }
  
  df.banff[is.na(df.banff)]=0
  
  df.banff <- apply(df.banff, 2, function(x) {x/sum(x)*100}) ## % biopsies
  
  suppressWarnings({
    df.banff <- reshape2::melt(df.banff)
  })
  df.banff$Var1 <- as.factor(df.banff$Var1)
  df.banff <- df.banff[order(df.banff$Var1, decreasing=TRUE),]
  df.banff[is.na(df.banff)]=0
  ## percent of patients (with a given banff score)
  ggplot(data=df.banff) + ylab("% biopsies") + ggtitle(cluster) +
    geom_bar(aes(x=Var2, y=value, fill=forcats::fct_rev(Var1)), color="black", stat="identity", alpha=1, width=0.9) + 
    scale_fill_manual("score", values = c("steelblue4", "steelblue3", "lightblue", "lightcyan")) +
    theme(axis.title.x=element_blank())
}

# g1 <- plot_k(aa.1.df, "k=1", index=TRUE)
# g2 <- plot_k(aa.2.df, "k=2", index=FALSE)
# g3 <- plot_k(aa.3.df, "k=3", index=TRUE)
# g4 <- plot_k(aa.4.df, "k=4", index=TRUE)
# grid.arrange(g1, g2, g3, g4, ncol=2, nrow=2)

## Optimal cluster # validation: silhouette index
## https://towardsdatascience.com/silhouette-coefficient-validating-clustering-techniques-e976bb81d10c


