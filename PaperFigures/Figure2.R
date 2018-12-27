  source('./Conf.R')
  source("./Main.R")
  
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  
  
  ######################################################################
  ##### 2C -- Association between folding energy and median PTR across tissues
  
  UTR5_FE_sliding <-readRDS(paste0(tempFileDir, "UTR5_FE_sliding.Rds"))
  
  #######Convert min.energy values in log2 scale in absolute value
  UTR5_FE_sliding[,-1] <- log2(abs(UTR5_FE_sliding[-1])+0.00001)
  
  RNA_Protein_CombinedUTR5Sliding <- merge(RNA_Protein_Combined,
                                           UTR5_FE_sliding, 
                                           by="EnsemblTranscriptID")
  
  allfeatures <- getAllFeatureColumns(RNA_Protein_CombinedUTR5Sliding)
  allfeatures$UTR5FE <- NULL
  selFeatures <- c(colnames(RNA_Protein_CombinedUTR5Sliding[,unlist(allfeatures)]))
  
  
  pVals <- c()
  betas <- c()
  for(elem in 1:(ncol(UTR5_FE_sliding)-1)){
    selFeaturesTemp <- c(paste0("Pos_",elem),selFeatures)
    lmRes <- lm(paste0("rowMedianPTR~", paste(selFeaturesTemp, collapse = "+")), 
                data=RNA_Protein_CombinedUTR5Sliding)
    PVALS <- summary(lmRes)$coefficients[,4]
    PVALS <- p.adjust(PVALS, method = "BH")
    pVals <- c(pVals,PVALS[2])
    betas <- c(betas, summary(lmRes)$coefficients[2,1])
  }
  
  UTR5SlidingDF <- data.frame(positionAsc=10^betas, 
                              pVals=pVals,
                              position=-75:75)
  UTR5SlidingDF$Significant <- F
  UTR5SlidingDF[UTR5SlidingDF$pVals< 0.1,"Significant"] <- T
  
  #### Fig 2C
  Figure2C <- ggplot(UTR5SlidingDF, aes(x=position, y=positionAsc, col=Significant)) +
    geom_point(size=0.5) + 
    theme(legend.position = "none")+
    xlab("Position relative to start codon") + 
    ylab("mRNA secondary structure effect \n on PTR ratio")+
    scale_y_continuous(trans = log2_trans(),breaks=seq(0.9,1.1,0.04), limits = c(0.9,1.1))+
    scale_x_continuous(breaks=seq(-75,75,10), limits = c(-75,75))+
    scale_color_manual(values=c("darkgrey", "red"))+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = 1, col="darkgrey")
  

  ##############################################################################
  ######### 2D uAUG effect
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                               responseStr=c("_PTR"))
  
  accTV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Motif_UTR5_AUG")
  
  Figure2D <- plotMotifEffectAcrossTissues(betas=accTV$betas, 
                                        tissueNames=TissueNamesPrint, 
                                        stdErrors=accTV$stdError,  
                                        ylabStr="out-frame uAUG effect \n on PTR ratio", 
                                        breaksSeq=seq(0.5,1,0.1), 
                                        labelsSeq=seq(0.5,1,0.1))

  
  ############################################################################
  #### Figure 2E
  
  accTV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Motif_UTR5_")
  
  pvals <- accTV$pvalues
  effectSizes <- 10^accTV$betas
  colnames(effectSizes) <- sapply(colnames(pvals), function(x){strsplit(x,"Motif_UTR5_")[[1]][2]})
  rownames(effectSizes) <- TissueNamesPrint
  effectSizes[pvals > 0.1] <- NA
  
  pheatmap(effectSizes,
           treeheight_row = 0,
           treeheight_col = 0,
           cluster_rows = F,
           cluster_cols = F,
           na_col = "grey",
           main="", legend = F, fontsize = 7)
  