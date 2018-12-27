  source('./Conf.R')
  source("./Main.R")

  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                               responseStr=c("_PTR"))
  
  length(which(RNA_Protein_Combined$PESTNoPotential>0))
  #############################################
  proteinFeatures <- getBetasOfFeature(allFeaturesFit,
                                       featureSubStr="ntermHydro|averagePISO|Motif_LinAA_|PESTNoPotential|Motif_AA_") 
  
  apply(10^proteinFeatures$betas,2, median )
  
  betas <- as.data.frame(10^proteinFeatures$betas)*100 -100
  
  betas$Tissues <- TissueNamesPrint
  pvals <- proteinFeatures$pvalues
  pvals$Tissues <- TissueNamesPrint
  
  meltedBetas <- melt(betas, id.vars="Tissues")
  meltedBetas$value <- meltedBetas$value -1
  
  
  meltedPvals <- melt(pvals, id.vars="Tissues")
  meltedPvals <- plyr::rename(meltedPvals, c("value" = "FDR"))
  meltedPvals$sign <- ""
  meltedPvals[meltedPvals$FDR < 0.1, "sign"] <- "*"
  
  aaMotifs <- c("KRR", "NS", "CG")
  linAAMotifs <- c("CLV_PCSK_FUR_1", "LIG_KEPE_1", "TRG_NLS_Bipartite_1", "TRG_NLS_MonoCore_2", 
                  "TRG_NLS_MonoExtC_3", "TRG_NLS_MonoExtN_4")
  
  meltedDF <- merge(meltedBetas, meltedPvals, by=c("Tissues","variable"))
  meltedDF$variable <- mapvalues(meltedDF$variable, from = levels(meltedDF$variable), to = c(aaMotifs,linAAMotifs, "PEST motif", "Protein isoelectric point", "Protein 5' end hydrophobicity"))
  
  
  featureColors <- c(rep("red", length(aaMotifs)), rep("navy", length(linAAMotifs)), "darkorange", "red4", "darkgreen")
  
  Fig4 <-  ggplot(data = meltedDF, aes(x=Tissues, y=variable,  fill=value)) +
    geom_tile() + 
    geom_text(aes(label = sign)) +
    scale_fill_gradient2(low = "blue", 
                        high = "darkred",
                        midpoint = 0,
                        mid="white",
                        limits=c(-70, 25))+
    ylab("")+xlab("")+
    theme(axis.text.x = element_text(family = fontFamily, size=7, angle = 90, hjust = 1),
          axis.text.y = element_text(family = fontFamily, size=7, colour = featureColors),
          axis.title.x=element_text(family = fontFamily, size=7), 
          axis.title.y=element_text(family = fontFamily, size=7),
          plot.title = element_text(family = fontFamily, size=7, hjust = 0.5))
