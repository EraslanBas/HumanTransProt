  source('./Conf.R')
  source("./Main.R")
  
  ################# 5A
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                               responseStr=c("_PTR"))
  
  accTV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Motif_UTR3_")
  
  pvals <- accTV$pvalues
  effectSizes <- 10^accTV$betas
  colnames(effectSizes) <- sapply(colnames(pvals), function(x){strsplit(x,"Motif_UTR3_")[[1]][2]})
  rownames(effectSizes) <- TissueNamesPrint
  effectSizes[pvals > 0.1] <- NA
  
  ## Figure 5A
  pheatmap(effectSizes,
           treeheight_row = 0,
           treeheight_col = 0,
           cluster_rows = F,
           cluster_cols = F,
           na_col = "grey",
           main="", legend = F, fontsize = 7)
  
  
  ################# 5C
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allFeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  
  X <- RNA_Protein_Combined[, unlist(allFeatures)]
  fixedEffectCols <- colnames(X)
  nFeature <- ncol(X)
  X_all <- X
  
  for(i in 1:28){
    X_all <- rbind(X_all, X)
  }
  
  Y <- RNA_Protein_Combined[,paste0(TissueNames,"_PTR")]
  Y <- melt(Y)
  Y <- plyr::rename(Y, c("variable"="tissuePTR", "value"="tissuePTRValue"))
  
  Z <- RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")]
  Z <- melt(Z)
  Z <- plyr::rename(Z,c("variable"="tissueRNA", "value"="tissueRNAValue"))
  
  
  X_all <- cbind(X_all, Y)
  X_all <- cbind(X_all, Z)
  X_all <- X_all[!is.na(X_all$tissuePTRValue),]
  
  res_lm <- compute10FCV_mixedEffectModel(X_all[,c(fixedEffectCols,"tissuePTR")], 
                                          y = X_all$tissuePTRValue, 
                                          rnaCol= X_all$tissueRNAValue, 
                                          fixedEffectCols = c(fixedEffectCols,"tissuePTR"), 
                                          randomEffectCol=NULL, 
                                          isMixedEffect=F, 
                                          flds = NULL)
  
  sumS <- sumStat(yValues=res_lm$yVal, yPredVlues=res_lm$yPred, includeCV=T)
  temp <- data.frame(yPred=10^res_lm$yPred, 
                     yVal=10^res_lm$yVal,
                     rnaVal=10^res_lm$rnaVal)
  
  Figure5C <- ggplotw(temp,aes(x=yPred,y=yVal)) + 
    geom_heatscatter(cor = T)+
    scale_x_continuous(trans = "log10", limits = c(10^3,10^7))+
    scale_y_continuous(trans = "log10", limits = c(1,10^9))+
    xlab("Predicted PTR ratio by features \n extracted from mRNA and protein sequences")+
    ylab("Observed PTR ratio ")+
    ggNT:::theme_Publication(base_family = 'Helvetica')+
    ourTheme+ theme(legend.position = c(0.7,0.1))
  

  ################# 5D
  
  Fig5D <- ggplotw(temp,aes(x=yPred,y=rnaVal)) + 
    geom_heatscatter(diag = F, cor=T)+
    ggNT:::theme_Publication(base_family = 'Helvetica')+
    scale_x_continuous(trans = "log10", limits = c(10^3,10^7))+
    scale_y_continuous(trans = "log10", limits = c(1,10^5))+
    xlab("Predicted PTR ratio by sequence features")+
    ylab("Observed mRNA levels")+
    ourTheme+
    theme(legend.position = c(0.4,0.1))
   
  
  ################# 5E
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allFeatures <- getAllFeatureColumnsWithProteinPTM(RNA_Protein_Combined)
  
  PTRexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined,
                                             isAAIncluded=T,
                                             responseSubstring="_PTR",
                                             ylabStr=" PTR ratio",
                                             allFeatures)
  
  
  ev <- 100*PTRexpVar$expVarsOfFeatures
  colnames(ev) <- featureCols[colnames(ev), "Feature"]
  

  ### Figure 5E
  pheatmap(ev, treeheight_row = 0, treeheight_col = 0, fontsize = 7)
  
  
  ################# 5F

  ptr<- plotResultFigures(expVarsOfFeatures=PTRexpVar$expVarsOfFeatures,
                          mainStr="Explained variance in protein to mRNA ratio (%)")
  ptr <- ptr$f
  colnames(ptr) <- c("PTR", "featureNames")
  ptr <- ptr[rownames(ptr) %ni% "ALL",]
  
  proHL <- readRDS(paste0(tempFileDir, "MathiesonProteinHLData.rds"))
  rownames(proHL) <- proHL$EnsemblGeneID
  proHL[,grep("Motif_LinAA_", colnames(proHL))] <- NULL
  
  kusterProteinHalfLife <- readRDS(paste0(tempFileDir, "/KusterProteinHalfLife.rds"))
  proHL <- merge(kusterProteinHalfLife, proHL, by="EnsemblGeneID", all.x=T, all.y=T)
  
  proHL <- proHL[,grep("_ProHL|ProteinHalfLife|EnsemblGeneID", colnames(proHL))]
  colnames(proHL) <- c("EnsemblGeneID", "BCells_ProHL", "Hepatocytes_ProHL", "Monocytes_ProHL")
  
  proHL <- merge(RNA_Protein_Combined, proHL, by="EnsemblGeneID")
  allFeatures <- getAllFeatureColumnsWithProteinPTM(proHL)
  
  ProHLexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined = proHL,
                                           isAAIncluded=T,
                                           responseSubstring=c("BCells_ProHL",
                                                               "Hepatocytes_ProHL",
                                                               "Monocytes_ProHL"),
                                           ylabStr=" protein half-life",
                                           allFeatures,
                                           isTissues=F)
  
  pro <- plotResultFigures(expVarsOfFeatures=ProHLexpVar$expVarsOfFeatures,
                           mainStr="Explained variance in protein half-life (%)")
  pro <- pro$f
  colnames(pro) <- c("ProteinHL", "featureNames")
  t <- as.data.frame(merge(ptr,pro, by="featureNames"))
  
  t$featureNames <- featureCols[t$featureNames, "Abreviation"]
  t[,c("PTR", "ProteinHL")] <- 100*t[,c("PTR", "ProteinHL")]
  cTest <-  cor.test(t$PTR, t$ProteinHL, method="spearman")
  p <- plotCorrelations(inputDF=t,
                        xStr="PTR",
                        yStr="ProteinHL",
                        xLab="Median explained variance in \n PTR ratio (%)",
                        yLab="Median explained variance in \n protein half-life (%)",
                        includeLabels = F)
  p <- p+ggtitle(paste(expression(rho), " = 0.59, P = 0.001"))
  Figure5F <- p + geom_text_repel(
    aes(label = featureNames),
    size = 2,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines")
  )
  
  
  
