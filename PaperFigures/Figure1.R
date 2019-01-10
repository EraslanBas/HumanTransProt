  source('./Conf.R')
  source("./Main.R")
  
  TissueAcro = c("Ad", "Ap", "Br", "Co", "Du", "En", "Es", "FT", "Fa", "GB", "He", "Ki", "Li", "Lu", "Ly", "Ov", "Pa", "Pl","Pr", "Re", "SG", "SI", "SM", "Sp", "St", "Te", "Th", "To", "UB")
  names(TissueAcro) = TissueNames
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  mRNA <- as.matrix(RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")])
  protein <- as.matrix(RNA_Protein_Combined[,paste0(TissueNames,"_iBAQ")])
  
  ## nested linear model
  ## P: matrix of proteins - or centered version of it
  ## P: matrix of mRNAs - or centered version of it
  nested_lm_res = function(R,P){
    as.data.table(t(
      sapply(
        1:ncol(P),
        function(j){
          y = P[,j]
          ## make sure to work with the same dataset for both fits for the test.
          any_na = is.na(y) | apply(is.na(R), 1, any) 
          y = P[!any_na,j]
          x = R[!any_na,j]
          x_rest = R[!any_na,-j]
          
          fit_single <- lm(y ~ x)
          fit_all <- lm(y ~ x + x_rest)
          p.value <- anova(fit_single, fit_all, test="Chisq")[["Pr(>Chi)"]][2]
          
          c( 
            r2_single_tissue = summary(lm(P[,j] ~ R[,j]))$adj.r.squared,
            r2_all_tissues = summary(lm(P[,j] ~R))$adj.r.squared,
            p.value = p.value
          )
        }
      )
    ))
  }
  
  ##### Figure 1B
  
  abs_vals_res = nested_lm_res(mRNA, protein)
  abs_vals_res[ ,Tissue:=TissueAcro[TissueNames]]
  
  abs_vals_res[which.min(r2_single_tissue)]
  abs_vals_res[which.max(r2_single_tissue)]
  
  
  abs_vals_res[which.min(r2_all_tissues)]
  abs_vals_res[which.max(r2_all_tissues)]
  
  Figure1B <- ggplot(abs_vals_res, aes(x=r2_single_tissue, y=r2_all_tissues))+
    geom_point(size=0.5, col="red")+
    geom_text_repel(aes(label=Tissue),
                    size=4)+
    theme(panel.spacing = unit(0, "lines"), 
          axis.text.x = element_text(size=12, hjust = 1, vjust=0.5), 
          axis.text.y = element_text(size=12, hjust = 1, vjust=0.5), 
          plot.title = element_text(size=12, hjust = 0.5), 
          legend.position = "None", strip.text = element_text(size=15))+
    geom_abline(intercept = 0, slope = 1)+
    xlim(c(0.2,0.5))+
    ylim(c(0.35,0.6))+
    xlab("Proportion of variance\n explained by mRNA of same tissue")+
    ylab("Proportion of variance\n explained by mRNA across tissues")

  
  ##### Figure 1C
  ## same for centered protein and mRNA levels
  
  mRNA_centered <- mRNA - rowMedians(mRNA, na.rm=TRUE)
  protein_centered <- protein - rowMedians(protein, na.rm=TRUE)
  rel_vals_res = nested_lm_res(mRNA_centered, protein_centered)
  rel_vals_res[ ,Tissue:=TissueAcro[TissueNames]]
  
  fig1C <- ggplot(rel_vals_res, aes(x=r2_single_tissue, y=r2_all_tissues))+
    geom_point(size=0.5, col="red")+
    geom_text_repel(aes(label=Tissue),
                    size=4)+
    theme(panel.spacing = unit(0, "lines"), 
          axis.text.x = element_text(size=12, hjust = 1, vjust=0.5), 
          axis.text.y = element_text(size=12, hjust = 1, vjust=0.5), 
          plot.title = element_text(size=12, hjust = 0.5), 
          legend.position = "None", strip.text = element_text(size=15))+
    geom_abline(intercept = 0, slope = 1)+
    xlim(c(0,0.6))+
    ylim(c(0,0.6))+
    xlab("Proportion of variance\n explained by mRNA of same tissue")+
    ylab("Proportion of variance\n explained by mRNA across tissues")
  
  
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
  
  
