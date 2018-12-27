  source('./Conf.R')
  source("./Main.R")

########### EV1 A, B
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[RNA_Protein_Combined$numNAPTR >5 ,]
  rownames(RNA_Protein_Combined) <- RNA_Protein_Combined$EnsemblGeneID
  
  tp <- t(RNA_Protein_Combined[RNA_Protein_Combined$numNAPTR==29,paste0(TissueNames,"_PTR")])
  rownames(tp) <- TissueNamesPrint
  
  pheatmap(tp,
           cluster_rows = F,
           treeheight_col = 0,
           show_colnames = F, fontsize = 7)

  tp <-apply(tp, 2, function(x){x-median(x, na.rm = T)})
  tp[tp>1] <- 1
  tp[tp< -1] <- -1
  
  pheatmap(tp,
           cluster_rows = F,
           treeheight_col = 0,
           show_colnames = F, fontsize = 7)

##########EV1 C,D,E,F
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[RNA_Protein_Combined$numNAPTR >5 ,]
  rownames(RNA_Protein_Combined) <- RNA_Protein_Combined$EnsemblGeneID
  
  PTR <- RNA_Protein_Combined[,paste0(TissueNames,"_PTR")]
  PTR <-t(apply(PTR, 1, function(x){x-median(x, na.rm = T)}))
  
  mRNA <- RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")]
  
  
  tp <- data.frame(medianMRNA=rowMedians(as.matrix(mRNA), na.rm = T),
                   ptrSD=apply(PTR,1,function(x){return(sd=sd(x, na.rm = T))}),
                   houseKeeping=RNA_Protein_Combined$housekeepingGene)
  p1 <- getDAVIDGO(foregroundGenes=rownames(PTR[tp$ptrSD < quantile(tp$ptrSD,c(0.1)),]), 
                   backgroundGenes=rownames(PTR),
                   fileNamePrefix="PTRLowVariance",
                   mainStr = "Genes with PTR standard deviation in the first decile")
  
  p1$plotList$GOTERM_BP_ALL
  p1$plotList$GOTERM_CC_ALL
  
  
  p2 <- getDAVIDGO(foregroundGenes=rownames(PTR[tp$ptrSD>quantile(tp$ptrSD,c(0.9)),]), 
                   backgroundGenes=rownames(PTR),
                   fileNamePrefix="PTRHighVariance",
                   mainStr = "Genes with PTR standard deviation in the last decile")
  
  p2$plotList$GOTERM_BP_ALL
  p2$plotList$GOTERM_CC_ALL
  
