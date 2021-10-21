
## predict probability of rejection on new biopsy
#wdir <- system("pwd", intern=TRUE)
#setwd(dirname(wdir))
#setwd(paste0(wdir, "/histomx/"))
source('../scripts/BHOT.R')

##-----------------------------------------------------
## load models
##-----------------------------------------------------
modelPath='../models/'

aa_model <- get(load(paste0(modelPath, 'aa_model.rda')))

normal.model <- get(load(paste0(modelPath, 'normal_model.RData')))
amr.model <- get(load(paste0(modelPath, 'amr_model.RData')))
tcmr.model <- get(load(paste0(modelPath, 'tcmr_model.RData')))
ati.model <- get(load(paste0(modelPath, 'ati_model.RData')))
ifta.model <- get(load(paste0(modelPath, 'ifta_model.RData')))
bkv.model <- get(load(paste0(modelPath, 'bkv_model.RData')))
g0_score.model <- get(load(paste0(modelPath, 'g0_score_model.RData')))
g1_score.model <- get(load(paste0(modelPath, 'g1_score_model.RData')))
ptc0_score.model <- get(load(paste0(modelPath, 'ptc0_score_model.RData')))
ptc1_score.model <- get(load(paste0(modelPath, 'ptc1_score_model.RData')))
i1_score.model <- get(load(paste0(modelPath, 'i1_score_model.RData')))
t1_score.model <- get(load(paste0(modelPath, 't1_score_model.RData')))
v0_score.model <- get(load(paste0(modelPath, 'v0_score_model.RData')))
cg0_score.model <- get(load(paste0(modelPath, 'cg0_score_model.RData')))
ci0_score.model <- get(load(paste0(modelPath, 'ci0_score_model.RData')))
ct0_score.model <- get(load(paste0(modelPath, 'ct0_score_model.RData')))
cv0_score.model <- get(load(paste0(modelPath, 'cv0_score_model.RData')))
c4d0_score.model <- get(load(paste0(modelPath, 'c4d0_score_model.RData')))

dx_normal <- c("Normal or minimal changes", "Pristine")
dx_amr <- c("Active AMR", "Chronic (+/- active) AMR")
dx_tcmr <- c("Acute TCMR", "Chronic active TCMR")
dx_ati <- c("Acute tubular injury")
dx_ifta <- c("Isolated IFTA") #ie. IFTA w/ no rejection
dx_viral <- c("NBKv")

## refset diagnosis and molecular score files
mscores_all_file <- '../static/refset_molecular_scores_all.txt'
mscores_aa_file <- '../static/refset_molecular_scores_aa.txt'

## import refSet molecular scores and Dx
mscores_ref <- read.table(mscores_all_file, header=TRUE, sep='\t')
rownames(mscores_ref) <- mscores_ref$ID
colnames(mscores_ref)[colnames(mscores_ref)=="amr"]<-"AMR"
colnames(mscores_ref)[colnames(mscores_ref)=="tcmr"]<-"TCMR"
colnames(mscores_ref)[colnames(mscores_ref)=="ati"]<-"ATI"
colnames(mscores_ref)[colnames(mscores_ref)=="ifta"]<-"IFTA"

##-------------
## reference set RCC files (for data normalization)
refRCCpath="../refRCCs/"
refRCC <- list.files(refRCCpath, pattern=".RCC", full.names=TRUE, recursive=TRUE)

## for testing
#newRCC="../test_files/test.RCC" 

##----------------------------------------
## single function for sample report
BHOTpred <- function(newRCC, outPath) {
  
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
  
  #cat("
  #     O       o O       o O       o
  #     | O   o | | O   o | | O   o |
  #     | | O | | | | O | | | | O | |
  #     | o   O | | o   O | | o   O |
  #     o       O o       O o       O\n")
  
  ## if no output directory provided, output files to RCC file directory
  if (exists('outPath')) {
    outPath=outPath
    cat("Output path:\n", outPath, "\n")
  } else {
    outPath=dirname(newRCC)
    cat("Output path:\n", outPath, "\n")
  }
  
  ## get new sample ID and verify BHOT sequencing panel
  new.ns.data <- parseRCC(newRCC)
  newID <- colnames(new.ns.data$attributes)[2]
  newRLF <- new.ns.data$attributes[new.ns.data$attributes$variable=="GeneRLF",newID]
  
  if (newRLF != "NS_Hs_Transplant_v1.0") {
    stop("The Gene RLF must be the same as the reference set: nCounter Human Organ Transplant Panel (NS_Hs_Transplant_v1.0)")
  } else {
    cat(">>Outputting results for sample", newID, "\n")
  }
  
  ## create sample directory in outPath
  dir.create(paste0(outPath, newID), recursive=TRUE, showWarnings=FALSE)
  newOut=paste0(outPath, newID)
  
  ## combine refset and new RCC into single count table
  RCCfiles <- c(refRCC, newRCC)
  
  ns.data <- parseRCC(RCCfiles) ##takes awhile
  
  countTable <- ns.data$counts
  rownames(countTable) <- countTable$Name
  
  ## exclude ref samples not used in classifiers
  countTable <- countTable[,colnames(countTable) %in% c("CodeClass", "Name", "Accession", rownames(mscores_ref), newID)]
  
  # exprs.all <- lapply(rccFiles, function(file) {
  #     countFile <- data.frame(nanostringr::parse_counts(file))
  #     return(countFile)
  # })
  #countTable <- Reduce(function(x,y) merge(x, y, all=FALSE, by=c("CodeClass", "Name", "Accession")), exprs.all)
  
  ## TODO: replace with manual normalization functions (see BHOT.R)
  new.ns.norm <- NanoStringNorm::NanoStringNorm(
            x = countTable,
            CodeCount = 'geo.mean', #pos controls
            Background = 'none', #mean.2sd
            SampleContent = 'housekeeping.geo.mean', #housekeeping genes
            round.values = FALSE, 
            take.log = TRUE,
            return.matrix.of.endogenous.probes=TRUE)
  
  ## new sample(s) normalized with refSet
  new.ns.norm <- data.frame(counts=new.ns.norm[,newID])
  colnames(new.ns.norm) <- newID
  
  ## TODO: prevent gene delimiter modification in model output
  rownames(new.ns.norm) <- gsub("-", ".", rownames(new.ns.norm))
  rownames(new.ns.norm) <- gsub("/", ".", rownames(new.ns.norm))
  rownames(new.ns.norm) <- gsub(" ", ".", rownames(new.ns.norm))
  #dim(new.ns.norm)
  
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
  ## predict ATI prob for new sample(s)
  #ati.genes <- rownames(ati.model$scaling) #LDA
  ati.genes <- names(ati.model$coefficients)[-1] #LR
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% ati.genes]
  
  #new.ati.pred <- predict(ati.model, newdata=new.ge, type="prob")$posterior #LDA
  new.ati.pred <- predict(ati.model, newdata=data.frame(t(new.ge)), type="response")
  
  new.ati.pred <- data.frame(new.ati.pred)
  #new.ati.pred$ID <- rownames(new.ati.pred)
  new.ati.pred$ID <- newID
  #new.ati.pred$other <- NULL #LDA
  colnames(new.ati.pred)<-c("ati", "ID")
  
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
  ## predict BKV prob for new sample(s) 
  # #bkv.genes <- rownames(bkv.model$scaling) #LDA
  # bkv.genes <- names(bkv.model$coefficients)[-1] #LR
  # 
  # new.ge <- new.ns.norm
  # new.ge <- t(new.ge)
  # new.ge <- new.ge[,colnames(new.ge) %in% bkv.genes]
  # 
  # #new.bkv.pred <- predict(bkv.model, newdata=new.ge, type="prob")$posterior #LDA
  # new.bkv.pred <- predict(bkv.model, newdata=data.frame(t(new.ge)), type="response") #LR
  # 
  # new.bkv.pred <- data.frame(new.bkv.pred)
  # new.bkv.pred$ID <- newID
  # #new.bkv.pred$other <- NULL #LDA
  # colnames(new.bkv.pred)<-c("bkv", "ID")
  
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
  ## predict prob g>1 for new sample(s)
  #g1.genes <- rownames(g1_score.model$scaling) #LDA
  g1.genes <- names(g1_score.model$coefficients)[-1] #LR
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% g1.genes]
  
  #new.pred <- predict(g1_score.model, newdata=new.ge, type="prob")$posterior #LDA
  new.pred <- predict(g1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  #new.pred$ID <- rownames(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.g1.pred <- new.pred
  colnames(new.g1.pred)<-c("g1_score", "ID")
  
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
  ## predict prob ptc>1 for new sample(s)
  #ptc1.genes <- rownames(ptc1_score.model$scaling) #LDA
  ptc1.genes <- names(ptc1_score.model$coefficients)[-1] #LR
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% ptc1.genes]
  
  #new.pred <- predict(ptc1_score.model, newdata=new.ge, type="prob")$posterior #LDA
  new.pred <- predict(ptc1_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  #new.pred$ID <- rownames(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.ptc1.pred <- new.pred
  colnames(new.ptc1.pred) <- c("ptc1_score", "ID")
  
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
  ## predict prob ci>0 for new sample(s)
  #ci0.genes <- rownames(ci0_score.model$scaling) #LDA
  ci0.genes <- names(ci0_score.model$coefficients)[-1] #LR 
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% ci0.genes]
  
  #new.pred <- predict(ci0_score.model, newdata=new.ge, type="prob")$posterior #LDA
  new.pred <- predict(ci0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.ci0.pred <- new.pred
  colnames(new.ci0.pred) <- c("ci0_score", "ID")
  
  ##--------------------------
  ## predict0 prob ct0>0 for new sample(s)
  #ct0.genes <- rownames(ct0_score.model$scaling) #LDA
  ct0.genes <- names(ct0_score.model$coefficients)[-1] #LR 
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% ct0.genes]
  
  #new.pred <- predict(ct0_score.model, newdata=new.ge, type="prob")$posterior #LDA
  new.pred <- predict(ct0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.ct0.pred <- new.pred
  colnames(new.ct0.pred) <- c("ct0_score", "ID")
  
  ##--------------------------
  ## predicv0 prob cv0>0 for new sample(s)
  #cv0.genes <- rownames(cv0_score.model$scaling) #LDA
  cv0.genes <- names(cv0_score.model$coefficients)[-1] #LR 
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% cv0.genes]
  
  #new.pred <- predict(cv0_score.model, newdata=new.ge, type="prob")$posterior #LDA
  new.pred <- predict(cv0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.cv0.pred <- new.pred
  colnames(new.cv0.pred) <- c("cv0_score", "ID")
  
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
  ## predict prob c4d>0 for new samples
  #c4d0.genes <- rownames(c4d0_score.model$scaling) #LDA
  c4d0.genes <- names(c4d0_score.model$coefficients)[-1] #LR
  
  new.ge <- new.ns.norm
  new.ge <- t(new.ge)
  new.ge <- new.ge[,colnames(new.ge) %in% c4d0.genes]
  
  #new.pred <- predict(c4d0_score.model, newdata=new.ge, type="prob")$posterior
  new.pred <- predict(c4d0_score.model, newdata=data.frame(t(new.ge)), type="response") #LR
  
  new.pred <- data.frame(new.pred)
  new.pred$ID <- newID
  new.pred$low <- NULL
  new.c4d0.pred <- new.pred
  colnames(new.c4d0.pred) <- c("c4d0_score", "ID")
  
  ##-------------------------------------
  ## score table
  options(scipen = 999)
  join_list <- list(new.normal.pred, new.amr.pred, new.tcmr.pred,
                    new.ati.pred, new.ifta.pred,
                    new.g0.pred, new.ptc0.pred, new.g1.pred, new.ptc1.pred,
                    new.i1.pred, new.t1.pred, new.v0.pred,
                    new.cg0.pred, new.cv0.pred, new.ci0.pred, new.ct0.pred)
  tab <- Reduce(function(...) merge(..., all=TRUE), join_list)
  rownames(tab) <- tab$ID
  
  colnames(tab)[colnames(tab)=="amr"]<-"AMR"
  colnames(tab)[colnames(tab)=="tcmr"]<-"TCMR"
  colnames(tab)[colnames(tab)=="ati"]<-"ATI"
  colnames(tab)[colnames(tab)=="ifta"]<-"IFTA"
  #colnames(tab)[colnames(tab)=="bkv"]<-"nBKV"
  
  mscores_new <- tab ## save scores before calculating %
  #tab <- data.frame(apply(tab[,-1], 2, function(x) {round(x, 3)*100})) #multiple samples
  tab <- data.frame(round(tab[-1], 3) *100) #single sample
  tab$ID <- rownames(tab)
  mscores_new_pct <- tab %>% dplyr::select(ID, everything())
  #write.table(mscores_new_pct, file=paste0(newOut, "/molecular_score_table_", newID, "_", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
  options(scipen = 0)
  
  ##--------------------------
  ## PCA of molecular scores
  ##--------------------------
  
  #table(mscores_ref$Dx)
  pca_scores <- c("normal", "AMR", "TCMR", "ATI", "IFTA") #pca input
  #pca_scores <- c("g0_score", "ptc0_score", "cg0_score", "i1_score", "t1_score")
  
  ## keep only Bx with all molecular scores
  mscores_pca <- mscores_ref[apply(mscores_ref[,pca_scores], 1, function(x) { all(!is.na(x)) }),]
  
  resPCA <- FactoMineR::PCA(mscores_pca[,pca_scores], scale.unit=FALSE, ncp=5, graph=FALSE)
  
  ## impute missing values
  #mc_scores_imputed <- missMDA::imputePCA(mc_scores, ncp=5, method="Regularized", scale=TRUE)
  #resPCA <- FactoMineR::PCA(mc_scores_imputed, scale.unit=FALSE, ncp=5, graph=FALSE)
  
  #factoextra::get_pca_var(resPCA)
  plot(resPCA, choix = "var") # variable graph
  #factoextra::fviz_screeplot(resPCA, ncp=5) #scree plot
  
  pca.df <- as.data.frame(resPCA$ind$coord)
  pca.df$ID <- rownames(pca.df)
  pca.df <- merge(pca.df, mscores_ref[,c("ID", "Dx")], by="ID", all.x=TRUE)
  
  ## simple Dx
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_amr, "AMR", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_tcmr, "TCMR", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_normal, "No specific Dx", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_ati, "ATI", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_ifta, "IFTA", pca.df$Dx)
  #pca.df$Dx <- ifelse(pca.df$Dx %in% dx_viral, "BKV", pca.df$Dx)
  
  #pca.df <- pca.df[pca.df$Dx %in% c("AMR", "TCMR", "No specific Dx"),]
  #my_cols=c("firebrick", "gray80", "steelblue")
  
  pca.df <- pca.df[pca.df$Dx %in% c("AMR", "TCMR", "No specific Dx", "ATI", "IFTA"),]
  my_cols=c("firebrick", "#009E73", "gray80", "black", "steelblue")
  
  ## colorblind friendly palette
  #my_cols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  #my_cols=c("firebrick", "gray80", "darkgray", "slategray", "turquoise", "paleturquoise", "forestgreen", "orange", "mediumpurple",
  #          "black", "steelblue", "springgreen", "ivory3", "salmon", "salmon", "dodgerblue", "darkviolet", "sienna", "blue3")
  #my_cols <- viridis::viridis(30)
  
  ## viral=PV, nBKV
  paste(levels(factor(pca.df$Dx)), my_cols)
  table(pca.df$Dx)
  
  ref_pca <- ggplot() +
    scale_fill_manual(values=my_cols) +
    geom_point(data=pca.df, aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) + 
    #geom_text_repel(data=pca.df, aes(Dim.1, Dim.2, label=ID), size=2, colour="gray", max.overlaps=10) +
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
          panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
  
  ## get outlier sample ID(s)
  #pca.df[pca.df$Dim.1 < -0.4 & pca.df$Dx != "TCMR",]
  #pca.df[pca.df$Dim.1 > 0.4 & pca.df$Dx != "AMR",]
  
  ##-----------------------
  ## output IQR and median scores by Dx for reference biopsies
  scores_keep <- c("AMR", "TCMR", "ATI", "IFTA", "normal",
                   "g0_score", "g1_score",
                   "ptc1_score", "ptc0_score", 
                   "i1_score", "t1_score", "v0_score",
                   "cg0_score", "ci0_score", "ct0_score", "cv0_score")
  
  #ref.tab <- merge(pca.df, mc_tab, by="ID") ## simplified Dx + PCs + molecular scores
  
  ref.tab <- mscores_ref ## complete Dx + molecular scores
  #ref.tab <- ref.tab[apply(ref.tab[,-1], 1, function(x) { all(!is.na(x)) }),] #exclude samples without all molecular scores
  
  ## simplify Dx
  ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_amr, "AMR", ref.tab$Dx)
  ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_tcmr, "TCMR", ref.tab$Dx)
  ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_normal, "No specific Dx", ref.tab$Dx)
  ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_ati, "ATI", ref.tab$Dx)
  ref.tab$Dx <- ifelse(ref.tab$Dx %in% dx_ifta, "IFTA", ref.tab$Dx)
  #table(ref.tab$Dx)
  
  ## IQR = Q3 - Q1
  ## Q0= min and Q5=max
  ref.amr.quantiles <- apply(ref.tab[ref.tab$Dx=="AMR",scores_keep], 2, quantile, na.rm=TRUE)
  ref.tcmr.quantiles <- apply(ref.tab[ref.tab$Dx=="TCMR",scores_keep], 2, quantile, na.rm=TRUE)
  ref.ati.quantiles <- apply(ref.tab[ref.tab$Dx=="ATI",scores_keep], 2, quantile, na.rm=TRUE)
  ref.ifta.quantiles <- apply(ref.tab[ref.tab$Dx=="IFTA",scores_keep], 2, quantile, na.rm=TRUE)
  ref.normal.quantiles <- apply(ref.tab[ref.tab$Dx=="No specific Dx",scores_keep], 2, quantile, na.rm=TRUE)
  #ref.bkv.quantiles <- apply(ref.tab[ref.tab$Dx=="nBKV",scores_keep], 2, quantile, na.rm=TRUE)
  
  ## score quartiles
  ref.amr.quantiles <- data.frame(apply(ref.amr.quantiles, 2, function(x) {round(x, 3)*100}))
  ref.amr.quantiles$Dx <- "AMR"
  ref.amr.quantiles$Q <- rownames(ref.amr.quantiles)
  
  ref.tcmr.quantiles <- data.frame(apply(ref.tcmr.quantiles, 2, function(x) {round(x, 3)*100}))
  ref.tcmr.quantiles$Dx <- "TCMR"
  ref.tcmr.quantiles$Q <- rownames(ref.tcmr.quantiles)
  
  ref.ati.quantiles <- data.frame(apply(ref.ati.quantiles, 2, function(x) {round(x, 3)*100}))
  ref.ati.quantiles$Dx <- "ATI"
  ref.ati.quantiles$Q <- rownames(ref.ati.quantiles)
  
  ref.ifta.quantiles <- data.frame(apply(ref.ifta.quantiles, 2, function(x) {round(x, 3)*100}))
  ref.ifta.quantiles$Dx <- "IFTA"
  ref.ifta.quantiles$Q <- rownames(ref.ifta.quantiles)
  
  ref.normal.quantiles <- data.frame(apply(ref.normal.quantiles, 2, function(x) {round(x, 3)*100}))
  ref.normal.quantiles$Dx <- "No specific Dx"
  ref.normal.quantiles$Q <- rownames(ref.normal.quantiles)
  
  ref.score.quantiles <- rbind(ref.amr.quantiles, ref.tcmr.quantiles, ref.ati.quantiles, ref.ifta.quantiles, ref.normal.quantiles)
  #write.table(ref.score.quantiles, file=paste0(newOut, "/refset_score_quantiles", Sys.Date(), ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
  
  ## resfset median scores (needed to reorder boxplots)
  ref.amr.median <- apply(ref.tab[ref.tab$Dx=="AMR",scores_keep], 2, median, na.rm=TRUE)
  ref.tcmr.median <- apply(ref.tab[ref.tab$Dx=="TCMR",scores_keep], 2, median, na.rm=TRUE)
  ref.ati.median <- apply(ref.tab[ref.tab$Dx=="ATI",scores_keep], 2, median, na.rm=TRUE)
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
  dat.m <- data.table::melt(dat)
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
  
  #ggplot(data=dat.m, aes(x=reorder(variable, value, median, order=TRUE), y=value)) + 
  boxplot_tcmr <- ggplot(data=dat.m, aes(x=variable, y=value)) + 
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
  #ggsave(paste0(newOut, "/tcmr_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_tcmr, device="pdf", width=8, height=6)
  
  ## ATI
  dat <- ref.tab[ref.tab$Dx=="ATI",]
  dat <- dat[,scores_keep]
  ati_order <- names(ref.ati.median[order(ref.ati.median, decreasing=F)])
  ## order by median scores (reorder not working)
  dat <- dat[,ati_order]
  colnames(dat) <- gsub("_score", "", colnames(dat))
  suppressWarnings({
  dat.m <- data.table::melt(dat)
  })
  dat.m$value <- dat.m$value * 100
  
  #ggplot(data=dat.m, aes(x=reorder(variable, value, median, order=TRUE), y=value)) + 
  boxplot_ati <- ggplot(data=dat.m, aes(x=variable, y=value)) + 
    ggtitle("ATI reference biopsies") +
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
  #ggsave(paste0(newOut, "/ati_reference_boxplot_", Sys.Date(), ".pdf"), plot=boxplot_tcmr, device="pdf", width=8, height=6)
  
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
  #all(colnames(mc_scores) == colnames(new.scoreTable[,scores_keep]))
  mscores_all <- rbind(mscores_ref[,pca_scores], mscores_new[,pca_scores])
  
  resPCA <- FactoMineR::PCA(mscores_all, scale.unit=FALSE, ncp=5, graph=FALSE)
  
  pca.df <- as.data.frame(resPCA$ind$coord)
  pca.df$ID <- rownames(pca.df)
  
  pca.df <- merge(pca.df, mscores_ref[,c("ID", "Dx")], by="ID", all.x=TRUE)
  pca.df$Dx <- ifelse(pca.df$ID %in% mscores_new$ID, "new", pca.df$Dx)
  
  #pca.df <- pca.df[!is.na(pca.df$Dx),]
  
  ## higlight new biopsy
  pca.df <- plyr::mutate(pca.df, ref=ifelse(pca.df$Dx!="new", "ref", "new"))
  
  ## simple Dx
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_amr, "AMR", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_tcmr, "TCMR", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_normal, "No specific Dx", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_ati, "ATI", pca.df$Dx)
  pca.df$Dx <- ifelse(pca.df$Dx %in% dx_ifta, "IFTA", pca.df$Dx)
  #pca.df$Dx <- ifelse(pca.df$Dx %in% dx_viral, "BKV", pca.df$Dx)
  
  #pca.df <- pca.df[pca.df$Dx %in% c("AMR", "TCMR", "No specific Dx", "new"),]
  pca.df <- pca.df[pca.df$Dx %in% c("AMR", "TCMR", "ATI", "IFTA", "No specific Dx", "new"),]
  
  table(pca.df$Dx)
  
  my_cols=c("firebrick", "#009E73", "gray80", "black", "steelblue")
  paste(levels(factor(pca.df$Dx)), my_cols)
  table(pca.df$Dx)
  
  pca_new_1_2 <- ggplot() +
    scale_fill_manual(values=my_cols) +
    #geom_point(data=pca.df, aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) + 
    geom_point(data=pca.df[pca.df$ref=="ref",], aes(Dim.1, Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) + 
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
          panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
  ggsave(paste0(newOut, "/pca_dx_pc1_pc2_", newID, "_", Sys.Date(), ".pdf"), plot=pca_new_1_2, device="pdf", width=7, height=6)
  
  pca_new_2_3 <- ggplot() +
    scale_fill_manual(values=my_cols) +
    geom_point(data=pca.df[pca.df$ref=="ref",], aes(x=Dim.3, y=Dim.2, fill=Dx), shape=21, color="gray", size=4, alpha=0.7) + 
    geom_point(data=pca.df[pca.df$ref=="new",], aes(x=Dim.3, y=Dim.2), shape=23, size=4, alpha=1, fill="orange") +
    geom_text_repel(data=pca.df[pca.df$ref=="new",], aes(x=Dim.3, y=Dim.2, label=ID), size=4, colour="orange") +
    xlab(paste("PC3 ", round(resPCA$eig[3,2]),"% of variance",sep="")) +
    ylab(paste("PC2 ", round(resPCA$eig[2,2]),"% of variance",sep="")) +
    theme(legend.position="top", legend.title=element_blank(),
          panel.grid.minor=element_line(colour="gray"), panel.grid.major=element_blank(),
          panel.background=element_blank(), 
          axis.line=element_line(colour="white"),
          axis.text=element_text(size=12), 
          axis.text.x=element_text(face="bold", size=12, angle=0),
          axis.title.x=element_text(face="bold", size=12, angle=0), 
          axis.title.y=element_text(face="bold", size=12),  
          panel.border=element_rect(colour="whitesmoke", fill=NA, size=1))
  #ggsave(paste0(newOut, "/pca_dx_pc2_pc3_", newID, "_", Sys.Date(), ".pdf"), plot=last_plot(), device="pdf", width=7, height=6)
  
  ##-----------------------------
  ## KNN: k-nearest neighbors
  ##-----------------------------
  k=50
  #knn_tab <- rbind(mscores_ref[,pca_scores], mscores_new[,pca_scores])
  knn_tab <- mscores_all
  knn_tab$ID <- rownames(knn_tab)
  knn_tab <- merge(knn_tab, pca.df[,c("ID", "Dx", "ref")], by="ID", all.x=TRUE)
  
  train_ref <- knn_tab[knn_tab$ref %in% "ref",c("normal", "AMR", "TCMR", "ATI", "IFTA", "ID", "Dx")]
  test_new <- knn_tab[knn_tab$ref %in% "new",c("normal", "AMR", "TCMR", "ATI", "IFTA")]
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
  length(nn_df$ID) == k
  nn_dx <- data.frame(table(knn_tab[knn_tab$ID %in% nn_df$ID,"Dx"]))
  colnames(nn_dx) <- c("Dx", "Total")
  nn_dx <- nn_dx[nn_dx$Dx!="new",]
  sum(nn_dx$Total) == k
  nn_dx$Percent <- round(nn_dx$Total/sum(nn_dx$Total) * 100, 2)
  
  nn_dx <- nn_dx[order(nn_dx$Total, decreasing=TRUE),]
  rownames(nn_dx) <- nn_dx$Dx
  nn_dx$Dx <- NULL
  #knitr::kable(t(nn_dx), row.names=TRUE, align='l', caption="", format="html") %>% kableExtra::kable_styling(position="left")
  
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
  ## archetypal analysis
  ##-----------------------------
  ## import cluster assignment for ref biopsies
  mscores_aa <- read.table(mscores_aa_file, header=TRUE, sep='\t')
  
  ## predict archetypes on unseen molecular scores
  mscores_new <- dplyr::rename(mscores_new, amr=AMR)
  mscores_new <- dplyr::rename(mscores_new, tcmr=TCMR)
  mscores_new <- dplyr::rename(mscores_new, ati=ATI)
  pred_aa <- data.frame(predict(aa_model, mscores_new[,colnames(aa_model$archetypes)]))
  rownames(pred_aa) <- newID
  mscores_new$aa_cluster <- gsub("X", "", colnames(pred_aa)[max.col(pred_aa)])
  mscores_new$Dx <- "new"

  mscores_aa_all <- rbind(mscores_aa, mscores_new[,colnames(mscores_aa)])
  
  pca.aa <- merge(pca.df, mscores_aa_all[,c("ID", "aa_cluster")], by="ID")
  pca.aa$aa_cluster <- as.factor(pca.aa$aa_cluster)
  
  pca.aa[pca.aa$ref=="new",]
  
  aa_cols=c("darkviolet", "black", "#009E73", "firebrick", "gray80", "steelblue")
  
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
  ggsave(paste0(newOut, "/pca_archetypes_pc1_2_", newID, "_", Sys.Date(), ".pdf"), plot=pca_aa, device="pdf", width=7, height=6)
  
  ggplot() +
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
  #ggsave(paste0(newOut, "/pca_archetypes_pc2_3_", newID, "_", Sys.Date(), ".pdf"), plot=last_plot(), device="pdf", width=7, height=6)
  
  ##---------------
  ## Dx by cluster
  ## banff score by cluster
  aa.df <- pca.aa[,c("ID", "aa_cluster", "Dx")]
  #aa.df <- aa.df[aa.df$Dx!="new",]
  
  aa.1.df <- droplevels(aa.df[aa.df$aa_cluster=="1",])
  aa.2.df <- droplevels(aa.df[aa.df$aa_cluster=="2",])
  aa.3.df <- droplevels(aa.df[aa.df$aa_cluster=="3",])
  aa.4.df <- droplevels(aa.df[aa.df$aa_cluster=="4",])
  
  ## Dx by cluster
  k1.df <- data.frame(table(aa.1.df$Dx))
  k1.df$k <- "1"
  colnames(k1.df) <- c("Dx", "total", "cluster")
  
  k2.df <- data.frame(k2=table(aa.2.df$Dx))
  k2.df$k <- "2"
  colnames(k2.df) <- c("Dx", "total", "cluster")
  
  k3.df <- data.frame(k3=table(aa.3.df$Dx))
  k3.df$k <- "3"
  colnames(k3.df) <- c("Dx", "total", "cluster")
  
  k4.df <- data.frame(k4=table(aa.4.df$Dx))
  k4.df$k <- "4"
  colnames(k4.df) <- c("Dx", "total", "cluster")
  
  ## Dx by cluster
  cluster_table <- dplyr::bind_rows(k1.df, k2.df, k3.df, k4.df)
  
  #write.table(cluster_table, file=paste0(newOut, "/aa_cluster_", newID, "_", Sys.Date(), ".txt"), row.names=F, quote=F, sep='\t')
  
  ##-------------------------------------------
  ## output files and plots for markdown report
  return(list(score_table=mscores_new_pct,
              ref_scores=ref.score.quantiles,
              knn=nn_dx,
              aa_cluster_table=cluster_table,
              aa_cluster_new=pred_aa, ##new sample cluster probs
              pca_molecular_scores=pca_new_1_2,
              pca_archetype=pca_aa))
  
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
  
  df.banff <- reshape2::melt(df.banff)
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


