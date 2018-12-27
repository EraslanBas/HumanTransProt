  source('./Conf.R')
  source("./Main.R")
  
  ##### Figure 1D  : range difference of mRNA and protein levels
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[RNA_Protein_Combined$numNAPTR >5 ,]
  rownames(RNA_Protein_Combined) <- RNA_Protein_Combined$EnsemblGeneID
  
  PTR <- RNA_Protein_Combined[,paste0(TissueNames,"_PTR")]
  PTR <-t(apply(PTR, 1, function(x){x-median(x, na.rm = T)}))
  
  mRNA <- RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")]
  
  
  tp <- data.frame(medianMRNA=rowMedians(as.matrix(mRNA), na.rm = T),
                   ptrSD=apply(PTR,1,function(x){return(sd=sd(x, na.rm = T))}),
                   houseKeeping=RNA_Protein_Combined$housekeepingGene)
  
  give.n <- function(x){
    return(c(y = median(x)+0.1, label = length(x))) 
  }
  
  p <- ggplot(tp,aes(x=houseKeeping, ptrSD, color=houseKeeping)) + 
    geom_boxplot(alpha=0.5) + 
    scale_color_manual(values = c('#999999','#E69F00')) +
    stat_summary(fun.data = give.n, geom = "text", fun.y = "median",position = position_dodge(width = 0.5))+
    theme(legend.position="None", legend.justification=c(0,1))+
    xlab("Housekeeping gene")+
    ylab("PTR ratio standard deviation \n across tissues (log10)")+
    labs(color='Housekeeping')  
  
  Figure1D <- add_pval(p,pairs = list(c(1,2)), heights = 1.5, log=T)
  

  ####### Figure 1E, F
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[RNA_Protein_Combined$numNAPTR >15 ,]
  rownames(RNA_Protein_Combined) <- RNA_Protein_Combined$EnsemblGeneID
  
  PTR <- RNA_Protein_Combined[,paste0(TissueNames,"_PTR")]
  mRNA <- RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")]
  rownames(PTR) <- RNA_Protein_Combined$EnsemblGeneID
  rownames(mRNA) <- RNA_Protein_Combined$EnsemblGeneID
  colnames(PTR) <- TissueNames
  colnames(mRNA) <- TissueNames
  
  PTR <- t(scale(t(PTR), center = T, scale = F))
  mRNA <- t(scale(t(mRNA), center = T, scale = F))
  
  tpData <- list(PTR=PTR, mRNA=mRNA)
  MOFAobject <- createMOFAobject(tpData)
  
  DataOptions <- getDefaultDataOptions()
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  TrainOptions <- getDefaultTrainOptions()
  MOFAobject <- prepareMOFA(
    MOFAobject,
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )
  
  MOFAobject <- runMOFA(MOFAobject, outfile="./DATA/MOFA_PTR_mRNA.hdf5")
  
  r2 <- calculateVarianceExplained(MOFAobject)
  
  ### Figure 1E
  barplot(r2$R2Total, ylim=c(0,0.65), col="darkgrey")

  
  ### Figure 1F
  pheatmap(r2$R2PerFactor, treeheight_row = 0, treeheight_col = 0, cluster_rows = F)
  
  
